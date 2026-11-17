"""
Reading the methods against each other.

Each method rests on different assumptions, so where they agree and where they
diverge says something about the microstructure. This module turns those
contrasts into explicit, checkable statements and a suggested headline value.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List

import numpy as np

from .correction import CorrectedPeak
from .methods import MethodResult


@dataclass
class Finding:
    level: str          # "critical", "warning", "info", "ok"
    topic: str
    text: str


@dataclass
class Interpretation:
    findings: List[Finding] = field(default_factory=list)
    headline: str = ""
    headline_size: float = np.nan
    headline_size_err: float = np.nan
    headline_strain: float = np.nan
    headline_strain_err: float = np.nan
    headline_method: str = ""
    headline_key: str = ""
    headline_size_sys: float = np.nan
    size_range: tuple = (np.nan, np.nan)
    strain_state: str = ""   # "resolved", "not_resolved", "inconsistent", "unknown"
    d_max: float = np.nan
    anisotropic: bool = False
    shape_determined: bool = True
    verdict: str = ""        # a few words: what kind of answer this is
    verdict_level: str = ""  # "ok", "info", "warning" or "critical"

    def add(self, level, topic, text):
        self.findings.append(Finding(level, topic, text))

    def set_verdict(self, text, level):
        self.verdict, self.verdict_level = text, level


def _t(val, err):
    if val is None or err is None or not np.isfinite(err) or err <= 0:
        return np.nan
    return val / err


def interpret(methods: Dict[str, MethodResult], peaks: List[CorrectedPeak],
              has_profile: bool, size_label: str = "L_vol", setup_critical=None) -> Interpretation:
    out = Interpretation()
    out = _interpret(out, methods, peaks, has_profile, size_label)
    if not has_profile and out.verdict_level != "critical" and out.verdict != "No peaks selected":
        out.set_verdict("Lower bound only: no instrument profile", "critical")
    for c in (setup_critical or []):
        out.findings.insert(0, Finding("critical", "Setup", c))
    if setup_critical:
        out.headline = ("UNRELIABLE — setup mismatch (see findings). No size or strain is reported until the "
                        "radiation / Kα2 setting matches the instrument profile.")
        out.headline_size = np.nan
        out.headline_key = ""
        out.set_verdict("Not reported: the settings do not match the instrument profile", "critical")
    return out


def _interpret(out, methods, peaks, has_profile, size_label):
    used = [p for p in peaks if p.use]
    valid = [p for p in used if not p.resolution_limited and p.beta_s > 0]
    n_rl = sum(p.resolution_limited for p in used)

    # ------------------------------------------------------------ data basis
    if not has_profile:
        out.add("critical", "Instrument",
                "No instrument profile: instrumental broadening is counted as sample broadening, so every "
                "size is a LOWER bound and strain is overestimated. Load a LaB₆ calibration.")
    if any("outside_calibrated_range" in p.flags for p in used):
        out.add("warning", "Instrument",
                "Some reflections lie outside the 2θ range of the standard; the instrument profile is extrapolated there.")
    if not used:
        out.headline = "No peaks selected."
        out.set_verdict("No peaks selected", "warning")
        return out
    if n_rl:
        lvl = "critical" if n_rl > len(used) / 2 else "warning"
        out.add(lvl, "Resolution",
                f"{n_rl} of {len(used)} reflections show no statistically resolved sample broadening; they give "
                "only lower bounds on size and are excluded from averages.")
    n_near = sum("near_resolution_limit" in p.flags and not p.resolution_limited for p in used)
    if n_near:
        out.add("warning", "Resolution",
                f"{n_near} reflection(s) are broadened by less than the resolution threshold relative to the "
                "instrument; their sizes are the least certain and most sensitive to the calibration.")
    if len(valid) >= 2:
        s = np.array([np.sin(p.theta) for p in valid])
        span = s.max() - s.min()
        if span < 0.25:
            out.add("warning", "Angular range",
                    f"sin θ spans only {span:.2f}: intercepts and slopes of the linear plots are poorly constrained.")
    n_ov = sum("overlapped" in p.fit.flags for p in used)
    if n_ov:
        out.add("info", "Overlap", f"{n_ov} reflection(s) were fitted as part of overlapping groups; their widths are less certain.")

    # ------------------------------------------------------------ all resolution-limited
    sch = methods.get("scherrer")
    if not valid:
        b = sch.size_bound_nm if sch else None
        out.headline = (f"Crystallite size is beyond the resolution of this instrument: {size_label} > {b:.0f} nm "
                        "(2σ lower bound). Mean sizes are not reported for such samples."
                        if b else "No resolvable sample broadening.")
        out.headline_method = "resolution bound"
        out.strain_state = "unknown"
        out.set_verdict("Lower bound only: beyond the instrument's resolution" if b
                        else "No resolvable sample broadening", "warning")
        return out

    # ------------------------------------------------------------ profile character
    lf = np.array([p.lorentz_fraction for p in valid if np.isfinite(p.lorentz_fraction)])
    n_bound = sum(p.fit.eta <= 0.005 or p.fit.eta >= 0.995 for p in valid)
    if len(lf) and n_bound <= len(valid) / 2:
        m = float(np.mean(lf))
        if m > 0.75:
            txt = f"Sample profiles are mostly Lorentzian (mean fraction {m:.2f}): typical of size broadening."
        elif m < 0.35:
            txt = f"Sample profiles are mostly Gaussian (mean fraction {m:.2f}): typical of strain / defect broadening."
        else:
            txt = f"Sample profiles are mixed (mean Lorentzian fraction {m:.2f}): size and strain both contribute."
        out.add("info", "Profile shape", txt)

    # ------------------------------------------------------------ is the profile SHAPE determined?
    # Single-line Voigt and Voigt W-H split each profile into Lorentzian (size) and Gaussian
    # (strain) parts using the fitted eta. If eta sits on a bound (0 or 1) for most reflections,
    # or is very uncertain, that split is not determined by the data and these two methods
    # must not decide anything.
    comp_methods = ("single_line", "voigt_wh")
    free = [p for p in valid if 0.005 < p.fit.eta < 0.995 and np.isfinite(p.fit.eta_err) and p.fit.eta_err < 0.25]
    if len(valid) >= 2 and len(free) < max(2, 0.5 * len(valid)):
        out.shape_determined = False
        n1 = sum(p.fit.eta >= 0.995 for p in valid)
        n0 = sum(p.fit.eta <= 0.005 for p in valid)
        tail = ""
        if n1 > len(valid) / 2:
            tail = (" Most sample peaks sit at η = 1: they are at least as heavy-tailed as a Lorentzian, which is "
                    "typical of a broad crystallite-size distribution (or of a mixture of sizes). Integral-breadth "
                    "sizes then assume the pseudo-Voigt shape and are approximate.")
        out.add("warning", "Profile shape",
                f"The Gaussian/Lorentzian split is not determined by the data (η at 1 for {n1}, at 0 for {n0}, "
                f"and poorly defined for {len(valid) - n1 - n0 - len(free)} of {len(valid)} reflections). "
                "Single-line Voigt and Voigt Williamson-Hall depend on that split, so they are not used for "
                "the strain decision or the headline." + tail)
        for key in comp_methods:
            r = methods.get(key)
            if r is not None:
                r.extra["shape_undetermined"] = True
                r.warnings.insert(0, "Not used for the interpretation: the sample profile shape (η) is not "
                                     "determined for most reflections.")
    use_comp = out.shape_determined

    # ------------------------------------------------------------ strain evidence
    # A strain counts as evidence only if it is statistically significant AND
    # material: its breadth 4 e tan(theta) must be >= 15 % of the sample breadth
    # at the highest angles. This guards against tiny but "significant" slopes
    # produced by small systematic errors in the profile model.
    sl = methods.get("single_line")
    top = sorted(valid, key=lambda p: p.theta)[-3:]

    def material(e):
        if e is None or not np.isfinite(e) or e <= 0:
            return False, 0.0
        f = float(np.mean([4 * e * np.tan(p.theta) / p.beta_s for p in top]))
        return f >= 0.15, f

    evidence = []
    for key in ("wh", "voigt_wh", "hw", "single_line"):
        r = methods.get(key)
        if key in comp_methods and not use_comp:
            continue
        if r and r.success and r.strain_significant is not None:
            sig = bool(r.strain_significant)
            if sig:
                mat, frac = material(r.strain)
                r.extra["strain_fraction_high_angle"] = frac
                if not mat:
                    out.add("info", "Strain",
                            f"{r.name}: strain {100 * r.strain:.3f} % is statistically significant but accounts for "
                            f"only {100 * frac:.0f} % of the high-angle breadth; treated as negligible.")
                    sig = False
            evidence.append((r.name, sig))
    if sch and sch.success and "r_sin" in sch.extra:
        rr, pp = sch.extra["r_sin"], sch.extra["p_sin"]
        trend = rr < -0.5 and pp < 0.05
        evidence.append(("Scherrer size vs angle", trend))
        out.add("info", "Size vs angle",
                f"Scherrer sizes vs sin θ: r = {rr:+.2f} (p = {pp:.2g}). "
                + ("Size falls with angle, as strain would cause." if trend else
                   "No significant fall with angle."))
    n_yes = sum(e for _, e in evidence)
    if evidence:
        names_yes = ", ".join(n for n, e in evidence if e) or "none"
        if n_yes == 0:
            out.strain_state = "not_resolved"
            out.add("ok", "Strain", "No method finds material strain. Size-only analysis (Scherrer) is appropriate.")
        elif n_yes / len(evidence) <= 0.25:
            out.strain_state = "not_resolved"
            out.add("ok", "Strain",
                    f"Only {names_yes} indicates strain; the other {len(evidence) - n_yes} tests do not. "
                    "Treated as not resolved (a lone small strain is typical of the approximations in that "
                    "method). Size-only analysis (Scherrer) is appropriate.")
        elif 2 * n_yes > len(evidence):
            out.strain_state = "resolved"
            out.add("info", "Strain",
                    f"Strain broadening is detected ({names_yes}). Scherrer sizes are therefore lower bounds; "
                    "use a size-strain method.")
        else:
            out.strain_state = "inconsistent"
            out.add("warning", "Strain",
                    f"Evidence for strain is mixed (material in: {names_yes}). The Scherrer size is quoted as the "
                    "headline but may be slightly low; treat strain values as indicative.")

    # ------------------------------------------------------------ anisotropy
    # Per-reflection sizes from Scherrer (if no strain) or single-line (if strain).
    # Anisotropy is claimed only if the spread is large (CV > 20 %) AND larger
    # than the per-reflection uncertainties allow (chi-square test, p < 0.01).
    src = sl if out.strain_state == "resolved" else sch
    if out.strain_state == "resolved" and not use_comp:
        src = None
        out.add("info", "Anisotropy",
                "Strain is present and the profile shape is not determined, so hkl-dependent size cannot be "
                "separated from strain here. The scatter of the points about the Halder-Wagner line is an upper "
                "limit on any anisotropy.")
    if src is not None and src.success:
        rows = [row for row in src.per_peak if row.get("size") and row.get("size_err")
                and np.isfinite(row["size_err"]) and row["size_err"] > 0]
        if len(rows) >= 3:
            from scipy.stats import chi2 as _chi2

            def _test(D, sD):
                w = 1 / sD**2
                mw = np.sum(w * D) / np.sum(w)
                c2 = float(np.sum(((D - mw) / sD) ** 2))
                dof = len(D) - 1
                return c2, dof, float(_chi2.sf(c2, dof)), float(100 * np.std(D, ddof=1) / np.mean(D)), mw

            D = np.array([row["size"] for row in rows])
            sD = np.array([row["size_err"] for row in rows])
            chi2, dof, p, cv, mw = _test(D, sD)
            src.extra["anisotropy_chi2_red"] = chi2 / dof
            src.extra["anisotropy_p"] = p
            # Robustness: the claim must survive dropping the single most deviant reflection,
            # so that one weak or overlapped peak cannot make a whole sample "anisotropic".
            single = None
            if cv > 20 and p < 0.01 and len(D) >= 4:
                k = int(np.argmax(np.abs(D - mw) / sD))
                keep = np.arange(len(D)) != k
                c2b, dofb, pb, cvb, _ = _test(D[keep], sD[keep])
                if not (cvb > 20 and pb < 0.01):
                    single = (rows[k], cvb, c2b / dofb)
            if cv > 20 and p < 0.01 and single is None:
                out.anisotropic = True
                why = ("without a strain trend" if out.strain_state != "resolved" else "after separating strain")
                out.add("warning", "Anisotropy",
                        f"Sizes differ by reflection (CV {cv:.0f}%, χ²/dof = {chi2 / dof:.1f}) {why}: crystallite "
                        "shape is anisotropic or there is hkl-dependent faulting. A single size is not "
                        "representative; report sizes per reflection (label them with hkl).")
            elif single is not None:
                row, cvb, c2b_red = single
                lab = f" ({row['hkl']})" if row.get("hkl") else ""
                out.add("warning", "Anisotropy",
                        f"Sizes differ by reflection (CV {cv:.0f}%, χ²/dof = {chi2 / dof:.1f}), but only because of "
                        f"the reflection at {row['two_theta']:.2f}°{lab} ({row['size']:.1f} nm); without it the "
                        f"CV is {cvb:.0f}% (χ²/dof = {c2b_red:.1f}), below the anisotropy criterion. Check that peak for weak intensity, "
                        "overlap or an impurity before concluding anything about crystallite shape; consider "
                        "unticking it.")
            elif cv > 20:
                out.add("info", "Anisotropy",
                        f"Per-reflection sizes scatter (CV {cv:.0f}%) but within their uncertainties "
                        f"(χ²/dof = {chi2 / dof:.1f}): noisy data rather than demonstrated anisotropy.")
            else:
                out.add("ok", "Anisotropy", f"Per-reflection sizes are consistent (CV {cv:.0f}%).")

    # ------------------------------------------------------------ Monshi slope
    mo = methods.get("monshi")
    if mo and mo.success and np.isfinite(mo.extra.get("slope_t", np.nan)):
        t = mo.extra["slope_t"]
        out.add("info" if abs(t) <= 2 else "warning", "Monshi slope",
                f"Monshi slope {mo.slope:.2f} ± {mo.slope_err:.2f} "
                + ("is consistent with 1 (size broadening)." if abs(t) <= 2 else
                   "is not 1: broadening is not size-only."))

    # ------------------------------------------------------------ linearity
    # Judge scatter about the line relative to the data level (R^2 is meaningless
    # for a flat, strain-free Williamson-Hall plot).
    for key in ("wh", "hw"):
        r = methods.get(key)
        if r and r.success and r.n >= 4 and r.plot.get("x"):
            x = np.array(r.plot["x"])
            y = np.array(r.plot["y"])
            yp = (r.slope * x + r.intercept) * r.plot.get("scale", 1.0)
            rel = float(np.sqrt(np.mean((y - yp) ** 2)) / np.mean(np.abs(y)))
            r.extra["relative_scatter"] = rel
            if rel > 0.15:
                out.add("warning", "Linearity",
                        f"{r.name}: points scatter by {100 * rel:.0f} % about the line; the isotropic model "
                        "does not describe this sample well.")

    # ------------------------------------------------------------ method agreement
    sizes = {r.name: r.size_nm for r in methods.values()
             if r.success and r.size_nm is not None and np.isfinite(r.size_nm) and r.size_nm > 0}
    if mo and mo.success and abs(mo.extra.get("slope_t", 0) or 0) > 2:
        sizes.pop(mo.name, None)   # Monshi size meaningless when slope != 1
    if not use_comp:
        for key in comp_methods:
            if methods.get(key) is not None:
                sizes.pop(methods[key].name, None)
    if sizes:
        v = np.array(list(sizes.values()))
        out.size_range = (float(v.min()), float(v.max()))
        spread = (v.max() - v.min()) / np.median(v)
        if len(v) >= 3:
            if spread < 0.15:
                out.add("ok", "Agreement", f"Methods agree on size within {100 * spread:.0f}% "
                        f"({v.min():.1f}–{v.max():.1f} nm).")
            else:
                out.add("info", "Agreement", f"Sizes span {v.min():.1f}–{v.max():.1f} nm across methods "
                        f"({100 * spread:.0f}% of the median). The spread reflects the different assumptions; "
                        "see the findings above for which assumptions hold.")

    if out.strain_state == "resolved" and sch and sch.success and sch.size_nm:
        sch.warnings.insert(0, "Lower bound: strain broadening is present and has been counted as size.")

    # ------------------------------------------------------------ measurable range
    # Largest size this instrument can resolve at mid-pattern angles (breadth ratio threshold).
    mids = [p for p in used if np.isfinite(p.beta_min) and p.beta_min > 0]
    d_max = None
    if mids and sch is not None:
        K_over = sch.extra.get("K_lam")
        d_max = float(np.median([K_over / (p.beta_min * np.cos(p.theta)) for p in mids])) if K_over else None
    if d_max and sch and sch.size_nm:
        out.d_max = d_max
        if sch.size_nm > 0.5 * d_max:
            out.add("warning", "Measurable range",
                    f"The Scherrer size ({sch.size_nm:.0f} nm) is more than half the largest size this instrument "
                    f"resolves (≈ {d_max:.0f} nm at the resolution threshold). Sizes this large are sensitive to the "
                    "calibration and should be quoted only within the lab's validated range.")
    b_rl = [row.get("size_bound") for row in (sch.per_peak if sch else []) if row.get("size_bound")]
    if b_rl and sch and sch.size_nm and max(b_rl) > sch.size_nm:
        out.add("warning", "Resolution",
                f"Unresolved reflections imply sizes above {max(b_rl):.0f} nm, larger than the mean of the resolved "
                f"ones ({sch.size_nm:.0f} nm): the mean is probably biased low by excluding them.")

    # ------------------------------------------------------------ beyond measurable range
    # If most reflections are within the resolution threshold (or unresolved),
    # a mean size is not trustworthy: report a conservative 2-sigma lower bound.
    n_weak = sum(p.resolution_limited or "near_resolution_limit" in p.flags for p in used)
    if sch is not None and n_weak > len(used) / 2:
        K_lam = sch.extra.get("K_lam")
        lbs = []
        for p in used:
            b = p.beta_upper if np.isfinite(p.beta_upper) and p.beta_upper > 0 else np.nan
            if K_lam and np.isfinite(b):
                lbs.append(K_lam / (b * np.cos(p.theta)))
        if lbs:
            lb = float(np.min(lbs))
            out.headline = (f"Beyond the reliable range of this instrument ({n_weak} of {len(used)} reflections are at or "
                            f"below the resolution threshold): {size_label} > {lb:.0f} nm (2σ lower bound). "
                            "Mean sizes are not reported for such samples.")
            out.headline_method = "resolution bound"
            out.strain_state = out.strain_state or "unknown"
            out.set_verdict("Lower bound only: beyond the instrument's reliable range", "warning")
            return out

    # ------------------------------------------------------------ headline
    hw = methods.get("hw")
    vwh = methods.get("voigt_wh") if use_comp else None
    if not use_comp:
        sl = None
    pick = None
    if out.anisotropic:
        src = sl if (out.strain_state == "resolved" and sl and sl.size_nm) else sch
        per = [row["size"] for row in src.per_peak if row.get("size")]
        out.headline = (f"Anisotropic: {size_label} = {min(per):.1f}–{max(per):.1f} nm depending on reflection "
                        f"({src.name}). Report per-reflection values.")
        out.headline_method = src.name
        out.headline_size = float(np.mean(per))
        out.set_verdict("Anisotropic: the size depends on direction", "warning")
        return out
    if out.strain_state in ("not_resolved", "inconsistent", "") and sch and sch.size_nm:
        pick = sch
    elif out.strain_state == "resolved" and not any(
            r is not None and r.success and r.size_nm and r.strain_significant for r in (vwh, hw)):
        # strain indicated (single-line / angle trend) but not quantified by a size-strain line
        pick = vwh if (vwh and vwh.success and vwh.size_nm) else (hw if hw and hw.size_nm else sch)
        out.headline_method = pick.name if pick else ""
        if pick is not None and pick.size_nm:
            out.headline_key = pick.key
            out.headline_size = pick.size_nm
            out.headline_size_err = pick.size_err if pick.size_err is not None else np.nan
            out.headline = (f"{size_label} ≈ {pick.size_nm:.1f} nm ({pick.name}); strain broadening is indicated "
                            "but not quantified reliably by the size-strain plots. Treat the size as approximate.")
            out.set_verdict("Approximate size: strain indicated but not quantified", "warning")
        else:
            out.headline = "Strain indicated but size and strain could not be separated; see per-method results."
            out.set_verdict("No single size: size and strain not separated", "warning")
        return out
    elif out.strain_state == "resolved":
        # Prefer Voigt W-H (separates the fitted Lorentzian and Gaussian components directly);
        # Halder-Wagner relies on a parabolic approximation to the Voigt breadth that inflates
        # small strains. Both are quoted so agreement can be seen.
        ok = lambda r: r is not None and r.success and r.size_nm and r.strain_significant  # noqa: E731
        pick = vwh if ok(vwh) else (hw if ok(hw) else (hw if hw and hw.size_nm else sch))
    else:
        pick = sch
    if pick is None or not pick.size_nm:
        out.headline = "No reliable single size; see per-method results."
        out.set_verdict("No single size", "warning")
        return out
    out.headline_method = pick.name
    out.headline_key = pick.key
    out.headline_size = pick.size_nm
    out.headline_size_err = pick.size_err if pick.size_err is not None else np.nan
    if pick.key == "scherrer":
        ul = [r.strain_upper for r in (hw, vwh) if r and r.success and r.strain_upper is not None]
        if out.strain_state == "inconsistent":
            st = "; strain evidence is mixed (see findings), so this size may be slightly low"
        else:
            ul = [u for u in ul if np.isfinite(u) and u > 0]
            st = (f"; strain not resolved (< {100 * max(ul):.3f} %)" if ul else "; strain negligible")
        n_txt = f"mean ± SD over {pick.n} reflections" if pick.n > 1 else "single reflection, ± propagated uncertainty"
        out.headline = (f"{size_label} = {pick.size_nm:.1f} ± {pick.size_err:.1f} nm "
                        f"(Scherrer, {n_txt}){st}.")
        if out.strain_state == "inconsistent":
            out.set_verdict("Size resolved; strain evidence mixed", "info")
        elif out.strain_state == "not_resolved":
            out.set_verdict("Size resolved; no strain resolved", "ok")
        else:
            out.set_verdict("Size resolved; strain not assessed", "info")
    else:
        out.headline_strain = pick.strain if pick.strain is not None else np.nan
        out.headline_strain_err = pick.strain_err if pick.strain_err is not None else np.nan
        err = f" ± {pick.size_err:.1f}" if pick.size_err else ""
        serr = f" ± {100 * pick.strain_err:.3f}" if pick.strain_err else ""
        out.headline = (f"{size_label} = {pick.size_nm:.1f}{err} nm, strain e = {100 * pick.strain:.3f}{serr} % "
                        f"({pick.name}, {pick.n} reflections).")
        other = hw if pick is vwh else None
        if other is not None and other.success and other.size_nm and other.strain is not None:
            out.headline += (f" Halder-Wagner: {other.size_nm:.1f} nm, {100 * other.strain:.3f} %.")
        if pick.strain_significant:
            out.set_verdict("Size and strain resolved", "ok")
        else:
            out.set_verdict("Size resolved; strain not significant", "info")
    if np.all(np.isfinite(out.size_range)) and out.size_range[1] > out.size_range[0]:
        out.headline += f" Across methods: {out.size_range[0]:.1f}–{out.size_range[1]:.1f} nm."
    if not has_profile:
        out.headline = "LOWER BOUND ONLY (no instrument profile). " + out.headline
        out.set_verdict("Lower bound only: no instrument profile", "critical")
    return out
