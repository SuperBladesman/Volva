"""
Instrument profile from a line-profile standard (e.g. NIST SRM 660c LaB6).

Model (default, "total"): the TOTAL FWHM of the standard peaks follows a
Caglioti curve and the Lorentzian fraction eta is constant, or varies linearly
with 2theta when that trend is statistically significant,

    FWHM^2 = U tan^2(theta) + V tan(theta) + W   (Caglioti, Paoletti & Ricci, 1958)
    eta    = e0 (+ e1 * 2theta)

Standard peaks are fitted as split (asymmetric) pseudo-Voigts, so axial-
divergence asymmetry at low angle does not distort their widths. Peaks whose
eta hits 0 or 1 are down-weighted, and eta is not extrapolated beyond the
calibrated range.

and the instrument's Gaussian (G) and Lorentzian (L) widths at any angle are
obtained by TCH decomposition of these SMOOTH curves. This is how Rietveld
programs describe a pseudo-Voigt instrument function. It is far more robust
than fitting G and L separately: the total width of each standard peak is well
determined, but its split into G and L is not (especially for sharp,
under-sampled peaks), and separate fits turn that scatter into zig-zags.

Legacy model ("components", Völva 4.1.0 profiles, still loadable): G^2 Caglioti and
L = X tan(theta) + Y / cos(theta), fitted to per-peak TCH components.

Profiles are saved as JSON so a lab can keep one calibrated profile per
instrument configuration.
"""

from __future__ import annotations

import datetime as _dt
import json
from dataclasses import asdict, dataclass, field
from itertools import product
from typing import List, Optional

import numpy as np

from . import __version__
from . import constants as C
from .background import auto_half_width, local_noise, noise_sigma, snip
from .io import Pattern
from .peaks import PeakFit, detect_peaks, fit_all
from .profiles import decompose, voigt_integral_breadth


@dataclass
class InstrumentProfile:
    name: str = "unnamed profile"
    created: str = ""
    volva_version: str = __version__
    wavelength: dict = field(default_factory=dict)
    doublet: bool = True
    standard: str = C.DEFAULT_STANDARD
    source_file: str = ""
    source_sha256: str = ""
    note: str = ""
    two_theta_range: tuple = (0.0, 180.0)
    U: float = 0.0
    V: float = 0.0
    W: float = 0.0
    X: float = 0.0
    Y: float = 0.0
    profile_model: str = "components"  # "total" (4.1.1+) or "components" (legacy)
    FU: float = 0.0                   # total-FWHM Caglioti (deg^2)
    FV: float = 0.0
    FW: float = 0.0
    f_model: str = "caglioti"         # or "poly2" (quadratic in 2theta) if Caglioti is not positive
    f_poly: list = field(default_factory=list)
    eta_coef: list = field(default_factory=list)   # eta = sum c_k (2theta)^k
    rms_F: float = 0.0
    rms_eta: float = 0.0
    r2_F: float = 0.0
    r2_eta: float = 0.0
    g_model: str = "caglioti"         # or "interp"
    l_model: str = "tan_sec"          # or "interp"
    rms_G: float = 0.0                # absolute RMS residual, degrees
    rms_L: float = 0.0
    r2_G: float = 0.0
    r2_L: float = 0.0
    peaks: list = field(default_factory=list)
    position_check: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)

    # ------------------------------------------------------------------
    def _theta(self, tt):
        return np.deg2rad(np.asarray(tt, float) / 2.0)

    # ---- total model -----------------------------------------------------
    def fwhm_at(self, tt):
        tt = np.asarray(tt, float)
        if self.f_model == "poly2" and self.f_poly:
            f = np.polyval(self.f_poly, tt)
            return np.clip(f, 1e-5, None)
        t = np.tan(self._theta(tt))
        return np.sqrt(np.clip(self.FU * t**2 + self.FV * t + self.FW, 1e-10, None))

    def eta_at(self, tt):
        # no extrapolation of the eta trend: held at its end values outside the calibrated range
        tt = np.clip(np.asarray(tt, float), *self.two_theta_range)
        c = self.eta_coef or [0.5]
        return np.clip(np.polyval(list(reversed(c)), tt), 0.0, 1.0)

    def _GL_total(self, tt, dF=0.0, dE=0.0):
        F = np.clip(self.fwhm_at(tt) + dF, 1e-6, None)
        E = np.clip(self.eta_at(tt) + dE, 0.0, 1.0)
        return decompose(F, E)

    # ---- public ------------------------------------------------------------
    def G_at(self, tt):
        tt = np.asarray(tt, float)
        if self.profile_model == "total":
            return self._GL_total(tt)[0]
        if self.g_model == "caglioti":
            t = np.tan(self._theta(tt))
            g2 = self.U * t**2 + self.V * t + self.W
            return np.sqrt(np.clip(g2, 0.0, None))
        xs = np.array([p["two_theta"] for p in self.peaks])
        ys = np.array([p["G"] for p in self.peaks])
        return np.interp(tt, xs, ys)

    def L_at(self, tt):
        tt = np.asarray(tt, float)
        if self.profile_model == "total":
            return self._GL_total(tt)[1]
        if self.l_model == "tan_sec":
            th = self._theta(tt)
            return np.clip(self.X * np.tan(th) + self.Y / np.cos(th), 0.0, None)
        xs = np.array([p["two_theta"] for p in self.peaks])
        ys = np.array([p["L"] for p in self.peaks])
        return np.interp(tt, xs, ys)

    def _total_sigmas(self, tt):
        # uncertainty of the fitted curves (standard error), not the scatter of single peaks
        n = max(sum(1 for p in self.peaks if p.get("used", True)), 3)
        sF = self.rms_F * np.sqrt(3.0 / n)
        sE = self.rms_eta * np.sqrt(2.0 / n)
        G0, L0 = self._GL_total(tt)
        Gf, Lf = self._GL_total(tt, dF=sF)
        Ge, Le = self._GL_total(tt, dE=sE)
        sG = np.sqrt((Gf - G0) ** 2 + (Ge - G0) ** 2)
        sL = np.sqrt((Lf - L0) ** 2 + (Le - L0) ** 2)
        return sG, sL

    def sigma_G(self, tt):
        g = self.G_at(tt)
        if self.profile_model == "total":
            return np.sqrt(self._total_sigmas(tt)[0] ** 2 + (C.INST_REL_UNC_FLOOR * g) ** 2)
        return np.sqrt(self.rms_G**2 + (C.INST_REL_UNC_FLOOR * g) ** 2)

    def sigma_L(self, tt):
        l_ = self.L_at(tt)
        if self.profile_model == "total":
            return np.sqrt(self._total_sigmas(tt)[1] ** 2 + (C.INST_REL_UNC_FLOOR * l_) ** 2)
        return np.sqrt(self.rms_L**2 + (C.INST_REL_UNC_FLOOR * l_) ** 2)

    def total_fwhm_at(self, tt):
        """Instrument FWHM (degrees) for either model."""
        if self.profile_model == "total":
            return self.fwhm_at(tt)
        from .profiles import tch_fwhm
        return tch_fwhm(self.G_at(tt), self.L_at(tt))

    def breadth_at(self, tt):
        return voigt_integral_breadth(self.G_at(tt), self.L_at(tt))

    def covers(self, tt, margin=2.0):
        lo, hi = self.two_theta_range
        return (lo - margin) <= tt <= (hi + margin)

    # ------------------------------------------------------------------
    def to_dict(self):
        d = asdict(self)
        d["two_theta_range"] = list(self.two_theta_range)
        return d

    def save(self, path: str):
        with open(path, "w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, indent=2, default=float)

    @classmethod
    def load(cls, path: str) -> "InstrumentProfile":
        with open(path, "r", encoding="utf-8") as f:
            d = json.load(f)
        d["two_theta_range"] = tuple(d.get("two_theta_range", (0, 180)))
        known = {k: v for k, v in d.items() if k in cls.__dataclass_fields__}
        prof = cls(**known)
        if "profile_model" not in d:
            prof.warnings.append("This profile uses the Völva 4.1.0 instrument model (separate G and L fits), "
                                 "which is less robust on real data. Recalibrate from the standard file.")
        return prof

    def summary(self) -> str:
        lines = [
            f"Instrument profile: {self.name}",
            f"  Standard: {self.standard}   file: {self.source_file}",
            f"  Created: {self.created}   Völva {self.volva_version}",
            f"  Wavelength Kα1 = {self.wavelength.get('ka1', float('nan')):.6f} nm, "
            f"doublet {'fitted' if self.doublet else 'not fitted'}",
            f"  Range: {self.two_theta_range[0]:.1f}–{self.two_theta_range[1]:.1f}° 2θ, "
            f"{len(self.peaks)} peaks",
        ]
        if self.profile_model == "total":
            if self.f_model == "caglioti":
                lines.append(f"  FWHM:  FWHM² = {self.FU:.3e} tan²θ + {self.FV:.3e} tanθ + {self.FW:.3e}  "
                             f"(R² = {self.r2_F:.3f}, RMS {self.rms_F:.4f}°)")
            else:
                lines.append(f"  FWHM:  quadratic in 2θ (Caglioti not positive)  (R² = {self.r2_F:.3f}, "
                             f"RMS {self.rms_F:.4f}°)")
            c = self.eta_coef
            lines.append(f"  η:     {c[0]:.3f} + {c[1]:.2e}·2θ  (RMS {self.rms_eta:.3f})" if len(c) > 1
                         else f"  η:     {c[0]:.3f} (constant)")
            lines.append("  G and L: TCH decomposition of the smooth FWHM and η curves")
        elif self.g_model == "caglioti":
            lines.append(f"  Gaussian:   G² = {self.U:.3e} tan²θ + {self.V:.3e} tanθ + {self.W:.3e}  "
                         f"(R² = {self.r2_G:.3f}, RMS {self.rms_G:.4f}°)")
        else:
            lines.append(f"  Gaussian:   interpolated between standard peaks (RMS {self.rms_G:.4f}°)")
        if self.profile_model != "total":
            if self.l_model == "tan_sec":
                lines.append(f"  Lorentzian: L = {self.X:.3e} tanθ + {self.Y:.3e} / cosθ  "
                             f"(R² = {self.r2_L:.3f}, RMS {self.rms_L:.4f}°)")
            else:
                lines.append(f"  Lorentzian: interpolated between standard peaks (RMS {self.rms_L:.4f}°)")
        pc = self.position_check
        if pc:
            lines.append(f"  Position check: {pc.get('message', '')}")
        for w in self.warnings:
            lines.append(f"  ⚠ {w}")
        if self.note:
            lines.append(f"  Note: {self.note}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
def _r2(y, yp):
    ss_tot = np.sum((y - y.mean()) ** 2)
    return float(1 - np.sum((y - yp) ** 2) / ss_tot) if ss_tot > 0 else 0.0


def _reflections(std_name, wl_nm, tt_max=160.0):
    info = C.STANDARDS.get(std_name)
    if not info:
        return None
    a = info["a_nm"]
    out = set()
    for h, k, l in product(range(0, 13), repeat=3):
        if (h, k, l) == (0, 0, 0) or not (h >= k >= l):
            continue
        if info["lattice"] == "cubic_Fd3m":
            allodd = h % 2 == k % 2 == l % 2 == 1
            alleven = h % 2 == k % 2 == l % 2 == 0 and (h + k + l) % 4 == 0
            if not (allodd or alleven):
                continue
        d = a / np.sqrt(h * h + k * k + l * l)
        s = wl_nm / (2 * d)
        if s < 1:
            tt = 2 * np.degrees(np.arcsin(s))
            if tt <= tt_max:
                out.add(round(tt, 4))
    return np.array(sorted(out))


def _position_check(fits_used, std_name, wl):
    ref = _reflections(std_name, wl["ka1"])
    if ref is None or len(fits_used) == 0:
        return {"message": "skipped"}
    obs = np.array([f.center for f in fits_used])
    d = np.array([obs_i - ref[np.argmin(np.abs(ref - obs_i))] for obs_i in obs])
    matched = np.abs(d) < 0.3
    frac = float(np.mean(matched))
    if frac < 0.6:
        return {"ok": False, "matched_fraction": frac,
                "message": (f"only {matched.sum()}/{len(obs)} peaks match {std_name} at "
                            f"λ = {wl['ka1']:.5f} nm. Wrong standard, wrong wavelength or large "
                            "zero error — check before using this profile.")}
    mean = float(np.mean(d[matched]))
    sd = float(np.std(d[matched]))
    ok = abs(mean) < 0.1
    msg = f"{matched.sum()}/{len(obs)} peaks indexed; mean Δ2θ = {mean:+.4f}° (SD {sd:.4f}°)"
    if not ok:
        msg += " — large offset: check zero error / specimen displacement"
    return {"ok": ok, "matched_fraction": frac, "mean_offset": mean, "sd_offset": sd, "message": msg}


def calibrate_from_standard(pattern: Pattern, wl: dict, doublet: bool = True,
                            standard: str = C.DEFAULT_STANDARD, name: str = "",
                            note: str = "", bg_half_width: Optional[float] = None,
                            min_r2: float = C.DEFAULT_STD_MIN_R2,
                            prominence_sigma: float = 10.0):
    """
    Build an InstrumentProfile from a standard pattern.

    Returns (profile, log_lines, fits, background).
    """
    log: List[str] = []
    tt, I = pattern.two_theta, pattern.intensity
    if bg_half_width is None:
        bg_half_width, note_bg = auto_half_width(tt, I)
        log.append(f"Background (SNIP): {note_bg}")
    bg = snip(tt, I, bg_half_width)
    ysub = I - bg
    sig_local = local_noise(I, bg)
    sig = float(np.median(sig_local))
    step = float(np.median(np.diff(tt)))
    cands = detect_peaks(tt, ysub, wl, doublet=doublet, prominence_sigma=prominence_sigma,
                         min_rel_height=0.005, noise=sig_local, log=log,
                         min_fwhm=lambda _tt: 1.5 * step)
    log.append(f"Standard: {len(cands)} candidate peaks")
    step = float(np.median(np.diff(tt)))
    fits = fit_all(tt, ysub, cands, wl, doublet=doublet, min_r2=min_r2, noise=sig_local, log=log,
                   min_fwhm=lambda _tt: 1.5 * step)

    used: List[PeakFit] = []
    for f in fits:
        why = f.message
        if f.accepted and "overlapped" in f.flags:
            f.accepted = False
            why = "overlapped with a neighbour"
        if f.accepted:
            used.append(f)
            log.append(f"  {f.center:8.3f}°: FWHM {f.fwhm:.4f}°  η {f.eta:.2f}  R² {f.r2:.4f}")
        else:
            log.append(f"  {f.center:8.3f}°: rejected ({why})")

    prof = InstrumentProfile(
        name=name or f"{pattern.name} ({standard})",
        created=_dt.datetime.now().isoformat(timespec="seconds"),
        wavelength=dict(wl), doublet=doublet, standard=standard,
        source_file=pattern.name, source_sha256=pattern.sha256, note=note)

    if len(used) < 3:
        prof.warnings.append(f"Only {len(used)} usable standard peaks; profile is unreliable.")
        log.append("✗ Too few standard peaks for a calibration.")
        if not used:
            return None, log, fits, bg

    ttp = np.array([f.center for f in used])
    F = np.array([f.fwhm for f in used])
    E = np.array([f.eta for f in used])
    Ee = np.array([f.eta_err if np.isfinite(f.eta_err) and f.eta_err > 0 else 0.05 for f in used])
    Ee = np.maximum(Ee, 0.02)                     # eta at a bound reports ~0 error
    G, L = decompose(F, E)
    theta = np.deg2rad(ttp / 2)
    t = np.tan(theta)
    grid = np.linspace(max(2.0, ttp.min() - 10), min(170.0, ttp.max() + 10), 200)
    tg = np.tan(np.deg2rad(grid / 2))

    # --- total FWHM: Caglioti with one round of outlier rejection ----------
    keep = np.ones(len(used), bool)
    cF = None
    for _ in range(3):
        if keep.sum() < 3:
            break
        A = np.column_stack([t**2, t, np.ones_like(t)])[keep]
        cF = np.linalg.lstsq(A, (F**2)[keep], rcond=None)[0]
        Fp = np.sqrt(np.clip(cF[0] * t**2 + cF[1] * t + cF[2], 1e-10, None))
        r = F - Fp
        mad = 1.4826 * np.median(np.abs(r[keep])) + 1e-5
        new_keep = np.abs(r) < 4 * mad
        if new_keep.sum() < 4 or np.array_equal(new_keep, keep):
            break
        for i in np.where(keep & ~new_keep)[0]:
            log.append(f"  {ttp[i]:8.3f}°: excluded from profile model (outlier, FWHM residual {r[i]:+.4f}°)")
        keep = new_keep

    prof.profile_model = "total"
    if cF is not None and np.all(cF[0] * tg**2 + cF[1] * tg + cF[2] > 0):
        prof.FU, prof.FV, prof.FW = map(float, cF)
        prof.f_model = "caglioti"
    else:
        deg = 2 if keep.sum() >= 4 else 1
        prof.f_poly = [float(c) for c in np.polyfit(ttp[keep], F[keep], deg)]
        prof.f_model = "poly2"
        prof.warnings.append("Caglioti curve for the total FWHM was not positive over the range; "
                             "a quadratic in 2θ is used instead.")

    # --- eta: weighted straight line in 2theta ------------------------------
    # A peak whose eta hit a bound (0 or 1) only tells us eta is "at least/at most" that
    # value, so it gets a large uncertainty rather than dragging the trend. The slope is
    # kept only if it is statistically significant; otherwise eta is a weighted constant.
    at_bound = (E >= 0.995) | (E <= 0.005)
    Ee = np.where(at_bound, np.maximum(Ee, C.ETA_BOUND_UNC), Ee)
    for i in np.where(at_bound & keep)[0]:
        log.append(f"  {ttp[i]:8.3f}°: η at bound ({E[i]:.2f}); down-weighted in the η trend")
    w = 1.0 / Ee[keep] ** 2
    const = float(np.sum(w * E[keep]) / np.sum(w))
    prof.eta_coef = [const]
    if keep.sum() >= 4:
        x = ttp[keep]
        A = np.column_stack([np.ones_like(x), x]) * np.sqrt(w)[:, None]
        coef, *_ = np.linalg.lstsq(A, E[keep] * np.sqrt(w), rcond=None)
        res_ = E[keep] - (coef[0] + coef[1] * x)
        dof = max(keep.sum() - 2, 1)
        s2 = max(float(np.sum(w * res_**2) / dof), 1.0)       # never smaller than the stated errors
        cov = np.linalg.inv(A.T @ A) * s2
        se1 = float(np.sqrt(cov[1, 1]))
        if abs(coef[1]) > 2 * se1:
            prof.eta_coef = [float(coef[0]), float(coef[1])]
            log.append(f"  η trend: slope {coef[1]:.2e} ± {se1:.1e} per degree (significant)")
        else:
            log.append(f"  η trend: slope {coef[1]:.2e} ± {se1:.1e} per degree not significant; "
                       f"constant η = {const:.3f} used")

    Fm = prof.fwhm_at(ttp)
    Em = prof.eta_at(ttp)
    prof.rms_F = float(np.sqrt(np.mean((F[keep] - Fm[keep]) ** 2)))
    prof.rms_eta = float(np.sqrt(np.mean((E[keep] - Em[keep]) ** 2)))
    prof.r2_F = _r2(F[keep], Fm[keep])
    prof.r2_eta = _r2(E[keep], Em[keep])

    prof.two_theta_range = (float(ttp.min()), float(ttp.max()))
    prof.peaks = [dict(two_theta=float(ttp[i]), fwhm=float(F[i]), eta=float(E[i]),
                       eta_err=float(Ee[i]), skew=float(used[i].skew),
                       G=float(G[i]), L=float(L[i]), r2=float(used[i].r2), used=bool(keep[i]))
                  for i in range(len(used))]
    # for information: how well the per-peak components follow the smooth model
    Gm, Lm = prof.G_at(ttp), prof.L_at(ttp)
    prof.rms_G = float(np.sqrt(np.mean((G[keep] - Gm[keep]) ** 2)))
    prof.rms_L = float(np.sqrt(np.mean((L[keep] - Lm[keep]) ** 2)))
    prof.r2_G = _r2(G[keep], Gm[keep])
    prof.r2_L = _r2(L[keep], Lm[keep])
    inside = np.linspace(ttp.min(), ttp.max(), 100)
    if np.any(prof.G_at(inside) <= 1e-4) or np.any(prof.L_at(inside) <= 1e-4):
        prof.warnings.append("The instrument profile becomes purely Gaussian or purely Lorentzian somewhere in the "
                             "range; check the η trend in the Instrument tab.")
    if prof.rms_F > 0.1 * float(np.median(F)):
        prof.warnings.append(f"Standard FWHMs scatter by {prof.rms_F:.4f}° about the Caglioti curve "
                             "(> 10 % of the typical width); check the standard fits.")
    step = pattern.step
    free = keep & ~at_bound
    eta_scatter = float(np.std(E[free])) if free.sum() >= 3 else 0.0
    if np.median(F) < 3 * step:
        prof.warnings.append(
            f"Standard peaks are only {np.median(F) / step:.1f} steps wide (FWHM {np.median(F):.3f}°, step "
            f"{step:.4f}°). Peak shape (η) is poorly defined; a finer step (≤ FWHM/5) is recommended "
            "for the calibration scan.")
    if at_bound[keep].sum() >= 2 or eta_scatter > 0.08:
        prof.warnings.append(
            f"Peak shape (η) of the standard is noisy (scatter {eta_scatter:.2f}; {int(at_bound[keep].sum())} "
            f"peak(s) at η = 0 or 1). A constant or smoothed η is used, but a calibration scan with a finer "
            f"step (≤ {np.median(F) / 5:.3f}°) and longer counting would define it better.")
    prof.position_check = _position_check(used, standard, wl)
    if prof.position_check.get("ok") is False:
        prof.warnings.append(prof.position_check["message"])
    log.append(prof.summary())
    return prof, log, fits, bg
