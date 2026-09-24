"""
Size and strain methods. All use integral breadths of the instrument-corrected
sample profile (radians, 2θ scale), the Kα1 wavelength and a shape constant K.

Strain convention (all methods): strain broadening beta_strain = 4 e tan(theta),
i.e. e is the Williamson-Hall "apparent strain". In reciprocal units
beta*_strain = 2 e d*.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np
from scipy import stats

from . import constants as C
from .correction import CorrectedPeak


@dataclass
class MethodResult:
    key: str
    name: str
    success: bool = False
    message: str = ""
    size_nm: Optional[float] = None
    size_err: Optional[float] = None
    size_bound_nm: Optional[float] = None        # lower bound when size not resolved
    size_spread: Optional[float] = None          # SD between peaks (per-peak methods)
    strain: Optional[float] = None               # fraction
    strain_err: Optional[float] = None
    strain_upper: Optional[float] = None         # upper limit when strain not resolved
    strain_significant: Optional[bool] = None
    r2: Optional[float] = None
    n: int = 0
    slope: Optional[float] = None
    slope_err: Optional[float] = None
    intercept: Optional[float] = None
    intercept_err: Optional[float] = None
    plot: dict = field(default_factory=dict)
    per_peak: List[dict] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    extra: dict = field(default_factory=dict)
    notes: str = ""
    kind: str = "linear"

    # convenience -------------------------------------------------------
    def size_text(self, digits=1):
        if self.size_nm is not None and np.isfinite(self.size_nm):
            e = f" ± {self.size_err:.{digits}f}" if self.size_err is not None and np.isfinite(self.size_err) else ""
            return f"{self.size_nm:.{digits}f}{e} nm"
        if self.size_bound_nm is not None and np.isfinite(self.size_bound_nm):
            return f"> {self.size_bound_nm:.0f} nm (not resolved)"
        return "—"

    def strain_text(self):
        if self.strain is None and self.strain_upper is None:
            return "—"
        if self.strain_significant is False and self.strain_upper is not None:
            if self.strain_upper <= 1e-7:
                return "none detected (no Gaussian broadening)"
            return f"< {100 * self.strain_upper:.3f} % (not resolved)"
        if self.strain is not None:
            e = f" ± {100 * self.strain_err:.3f}" if self.strain_err is not None and np.isfinite(self.strain_err) else ""
            return f"{100 * self.strain:.3f}{e} %"
        return "—"


# ---------------------------------------------------------------------------
# regression helpers
# ---------------------------------------------------------------------------
def ols(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    A = np.column_stack([x, np.ones(n)])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    yp = A @ coef
    ss_res = float(np.sum((y - yp) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    if n > 2:
        s2 = ss_res / (n - 2)
        cov = s2 * np.linalg.inv(A.T @ A)
        se = np.sqrt(np.diag(cov))
    else:
        se = np.array([np.nan, np.nan])
    return float(coef[0]), float(coef[1]), float(se[0]), float(se[1]), r2


def ols_origin(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    b = float(np.sum(x * y) / np.sum(x * x))
    res = y - b * x
    se = float(np.sqrt(np.sum(res**2) / (n - 1) / np.sum(x * x))) if n > 1 else np.nan
    return b, se


def _valid(peaks):
    return [p for p in peaks if p.use and not p.resolution_limited and p.beta_s > 0]


def _bound(p, K, lam):
    """Lower size bound for one unresolved reflection, from its 2-sigma breadth upper limit."""
    b = p.beta_upper if np.isfinite(p.beta_upper) and p.beta_upper > 0 else p.beta_min
    return float(K * lam / (b * np.cos(p.theta))) if np.isfinite(b) and b > 0 else np.nan


def _bound_from(peaks, K, lam):
    """Lower size bound from resolution-limited peaks (true in every measured direction)."""
    b = [_bound(p, K, lam) for p in peaks if p.use and p.resolution_limited]
    b = [x for x in b if np.isfinite(x)]
    return float(np.min(b)) if b else None


def _pct_err(samples):
    s = samples[np.isfinite(samples)]
    if len(s) < 10:
        return np.nan
    lo, hi = np.percentile(s, [16, 84])
    return float(0.5 * (hi - lo))


# ---------------------------------------------------------------------------
NOTES = {
    "scherrer": (
        "Assumes all sample broadening is due to crystallite size (no strain). With integral "
        "breadth and K = 1 each value is a volume-weighted column length along that reflection's "
        "normal. If strain is present the values are LOWER bounds and fall with angle. "
        "Spread between reflections reflects anisotropy, strain or faulting."),
    "single_line": (
        "Single-line Voigt method (de Keijser et al., J. Appl. Cryst. 15, 308, 1982): per "
        "reflection, the Lorentzian part of the sample profile is taken as size and the Gaussian "
        "part as strain. Needs good counting statistics; fragile when the sample broadening is "
        "small compared with the instrument."),
    "monshi": (
        "Monshi, Foroughi & Monshi (World J. Nano Sci. Eng. 2, 154, 2012): least-squares line "
        "through ln β vs ln(1/cos θ). For pure size broadening the slope is 1; a slope "
        "significantly above 1 indicates strain (or other angle-dependent broadening), and the "
        "intercept-derived size is then not meaningful."),
    "wh": (
        "Williamson-Hall, uniform deformation model (Williamson & Hall, Acta Metall. 1, 22, 1953). "
        "Assumes size and strain breadths add linearly, which is strictly true only if both are "
        "Lorentzian. For Voigt profiles (Lorentzian size, Gaussian strain) the combined breadth is "
        "smaller than the linear sum, so the fitted slope (strain) is typically too small and the "
        "size too large; compare with Halder-Wagner and Voigt W-H. Isotropic strain assumed."),
    "voigt_wh": (
        "Voigt Williamson-Hall (after Langford, NIST Spec. Publ. 846, 110, 1992): separate plots "
        "for the Lorentzian breadths (size from the intercept) and Gaussian breadths (strain "
        "from the slope through the origin). Uses the profile shape rather than assuming how the "
        "breadths combine."),
    "hw": (
        "Halder-Wagner (Acta Cryst. 20, 312, 1966): Voigt approximation β*² = β*·β*_L + β*_G², "
        "plotted as (β*/d*)² vs β*/d*². Down-weights high-angle reflections. The size-strain "
        "plot (SSP), written in its dimensionally consistent form, is algebraically identical "
        "(Pal, J. Alloys Compd. 1048, 185324, 2025) and is therefore not reported separately."),
}


def scherrer(peaks: List[CorrectedPeak], lam, K) -> MethodResult:
    r = MethodResult("scherrer", "Scherrer (integral breadth)", kind="per_peak", notes=NOTES["scherrer"])
    used = [p for p in peaks if p.use]
    valid = _valid(peaks)
    valid_ids = {id(p) for p in valid}
    for p in used:
        row = dict(two_theta=p.two_theta, hkl=p.hkl, resolution_limited=p.resolution_limited)
        if id(p) in valid_ids:
            D = K * lam / (p.beta_s * np.cos(p.theta))
            with np.errstate(divide="ignore", invalid="ignore"):
                Dmc = K * lam / (p.mc["beta_s"] * np.cos(p.theta))
            row.update(size=D, size_err=_pct_err(np.where(p.mc["beta_s"] > 0, Dmc, np.nan)))
        else:
            row.update(size=None, size_bound=_bound(p, K, lam))
        r.per_peak.append(row)
    r.n = len(valid)
    if not used:
        r.message = "No peaks selected."
        return r
    r.extra["K_lam"] = K * lam
    if not valid:
        r.success = True
        r.size_bound_nm = _bound_from(used, K, lam)
        r.message = "All reflections are resolution-limited."
        r.warnings.append("Crystallites are too large for this instrument: only a lower bound can be given.")
        return r
    D = np.array([row["size"] for row in r.per_peak if row.get("size") is not None])
    s = np.array([np.sin(p.theta) for p in valid])
    r.success = True
    r.extra["K_lam"] = K * lam
    r.size_nm = float(np.mean(D))
    r.size_spread = float(np.std(D, ddof=1)) if len(D) > 1 else np.nan
    # one reflection: fall back to its own propagated uncertainty
    r.size_err = r.size_spread if len(D) > 1 else float(r.per_peak[[i for i, row in enumerate(r.per_peak)
                                                                     if row.get("size") is not None][0]]["size_err"])
    r.extra["cv_percent"] = float(100 * r.size_spread / r.size_nm) if len(D) > 1 else np.nan
    r.extra["median"] = float(np.median(D))
    r.extra["min"], r.extra["max"] = float(D.min()), float(D.max())
    if len(D) >= 4:
        with np.errstate(all="ignore"):
            import warnings as _w
            with _w.catch_warnings():
                _w.simplefilter("ignore")
                rr, pp = stats.pearsonr(s, D)
        r.extra["r_sin"] = float(rr)
        r.extra["p_sin"] = float(pp)
    n_rl = sum(p.resolution_limited for p in used)
    if n_rl:
        r.warnings.append(f"{n_rl} resolution-limited reflection(s) excluded from the mean.")
    r.plot = dict(x=[p.two_theta for p in valid], y=list(D),
                  yerr=[row["size_err"] for row in r.per_peak if row.get("size") is not None],
                  xb=[row["two_theta"] for row in r.per_peak if row.get("size") is None],
                  yb=[row.get("size_bound") for row in r.per_peak if row.get("size") is None],
                  xlabel="2θ (°)", ylabel="Size (nm)")
    return r


def single_line(peaks, lam, K) -> MethodResult:
    r = MethodResult("single_line", "Single-line Voigt", kind="per_peak", notes=NOTES["single_line"])
    valid = _valid(peaks)
    if not valid:
        r.message = "No resolved reflections."
        r.size_bound_nm = _bound_from(peaks, K, lam)
        r.success = r.size_bound_nm is not None
        return r
    Ds, Derr, es, eerr, sig_L, sig_G = [], [], [], [], [], []
    for p in valid:
        c, t = np.cos(p.theta), np.tan(p.theta)
        L_ok = p.beta_Ls > 2 * p.beta_Ls_err
        G_ok = p.beta_Gs > 2 * p.beta_Gs_err
        row = dict(two_theta=p.two_theta, hkl=p.hkl, lorentz_fraction=p.lorentz_fraction)
        if L_ok:
            D = K * lam / (p.beta_Ls * c)
            Dmc = K * lam / (np.where(p.mc["beta_Ls"] > 0, p.mc["beta_Ls"], np.nan) * c)
            row.update(size=D, size_err=_pct_err(Dmc))
            Ds.append(D)
            Derr.append(row["size_err"])
        else:
            den = (p.beta_Ls + 2 * p.beta_Ls_err) * c
            row.update(size=None, size_bound=(K * lam / den) if den > 0 else None)
        e = p.beta_Gs / (4 * t)
        e_err = p.beta_Gs_err / (4 * t)
        row.update(strain=e, strain_err=e_err, strain_resolved=bool(G_ok))
        es.append(e)
        eerr.append(e_err)
        sig_L.append(L_ok)
        sig_G.append(G_ok)
        r.per_peak.append(row)
    r.n = len(valid)
    r.success = True
    if Ds:
        r.size_nm = float(np.mean(Ds))
        r.size_spread = float(np.std(Ds, ddof=1)) if len(Ds) > 1 else np.nan
        r.size_err = r.size_spread if len(Ds) > 1 else float(Derr[0])
        r.extra["cv_percent"] = float(100 * r.size_spread / r.size_nm) if len(Ds) > 1 else np.nan
        if len(Ds) < len(valid):
            r.warnings.append(f"Lorentzian (size) broadening not resolved for {len(valid) - len(Ds)} reflection(s).")
    else:
        bounds = [row["size_bound"] for row in r.per_peak if row.get("size_bound")]
        r.size_bound_nm = float(np.min(bounds)) if bounds else None
        r.warnings.append("No reflection shows resolvable Lorentzian (size) broadening.")
    es = np.array(es)
    r.strain = float(np.mean(es))
    r.strain_err = float(np.std(es, ddof=1)) if len(es) > 1 else float(eerr[0])
    # Isotropic strain should give similar e for every reflection. Count it as
    # resolved only if most reflections show Gaussian broadening, the mean is
    # > 2 standard errors, and the per-reflection values agree (CV < 50 %).
    se = r.strain_err / np.sqrt(len(es)) if len(es) > 1 else r.strain_err
    cv_e = r.strain_err / r.strain if r.strain > 0 else np.inf
    r.extra["strain_cv"] = float(cv_e)
    most = np.sum(sig_G) > len(sig_G) / 2
    r.strain_significant = bool(most and r.strain > 2 * se and cv_e < 0.5)
    if most and cv_e >= 0.5:
        r.warnings.append(f"Per-reflection strains are inconsistent (CV {100 * cv_e:.0f} %): the Gaussian "
                          "components are poorly determined, so single-line strain is not reliable here.")
    if not r.strain_significant:
        r.strain_upper = float(np.mean(es + 2 * np.array(eerr)))
    r.plot = dict(x=[row["two_theta"] for row in r.per_peak if row.get("size") is not None],
                  y=Ds, yerr=Derr, xs=[row["two_theta"] for row in r.per_peak],
                  ys=list(100 * es), yserr=list(100 * np.array(eerr)),
                  xlabel="2θ (°)", ylabel="Size (nm)")
    return r


def monshi(peaks, lam, K) -> MethodResult:
    r = MethodResult("monshi", "Monshi-Scherrer", notes=NOTES["monshi"])
    v = _valid(peaks)
    r.n = len(v)
    if len(v) < C.MIN_PEAKS_MONSHI:
        r.message = f"Needs ≥ {C.MIN_PEAKS_MONSHI} resolved reflections (have {len(v)})."
        return r
    x = np.array([np.log(1 / np.cos(p.theta)) for p in v])
    y = np.array([np.log(p.beta_s) for p in v])
    yerr = np.array([p.beta_s_err / p.beta_s for p in v])
    m, b, sm, sb, r2 = ols(x, y)
    r.success = True
    r.slope, r.intercept, r.slope_err, r.intercept_err, r.r2 = m, b, sm, sb, r2
    r.size_nm = float(K * lam / np.exp(b))
    r.size_err = float(r.size_nm * sb) if np.isfinite(sb) else None
    t = (m - 1) / sm if sm and np.isfinite(sm) and sm > 0 else np.nan
    r.extra["slope_t"] = float(t)
    if np.isfinite(t) and abs(t) > 2:
        r.warnings.append(
            f"Slope {m:.2f} ± {sm:.2f} differs significantly from 1: broadening is not pure size, "
            "so the Monshi size is not physically meaningful.")
    r.plot = dict(x=list(x), y=list(y), yerr=list(yerr), xlabel="ln(1/cos θ)", ylabel="ln β (β in rad)")
    return r


def williamson_hall(peaks, lam, K) -> MethodResult:
    r = MethodResult("wh", "Williamson-Hall (UDM)", notes=NOTES["wh"])
    v = _valid(peaks)
    r.n = len(v)
    if len(v) < C.MIN_PEAKS_LINEAR:
        r.message = f"Needs ≥ {C.MIN_PEAKS_LINEAR} resolved reflections (have {len(v)})."
        return r
    x = np.array([4 * np.sin(p.theta) for p in v])
    y = np.array([p.beta_s * np.cos(p.theta) for p in v])
    yerr = np.array([p.beta_s_err * np.cos(p.theta) for p in v])
    m, b, sm, sb, r2 = ols(x, y)
    r.success = True
    r.slope, r.intercept, r.slope_err, r.intercept_err, r.r2 = m, b, sm, sb, r2
    _size_from_intercept(r, b, sb, K * lam)
    _strain_from_slope(r, m, sm)
    r.plot = dict(x=list(x), y=list(y * 1e3), yerr=list(yerr * 1e3),
                  xlabel="4 sin θ", ylabel="β cos θ (mrad)", scale=1e3)
    return r


def voigt_wh(peaks, lam, K) -> MethodResult:
    r = MethodResult("voigt_wh", "Voigt Williamson-Hall", notes=NOTES["voigt_wh"])
    v = _valid(peaks)
    r.n = len(v)
    if len(v) < C.MIN_PEAKS_LINEAR:
        r.message = f"Needs ≥ {C.MIN_PEAKS_LINEAR} resolved reflections (have {len(v)})."
        return r
    x = np.array([4 * np.sin(p.theta) for p in v])
    yL = np.array([p.beta_Ls * np.cos(p.theta) for p in v])
    yG = np.array([p.beta_Gs * np.cos(p.theta) for p in v])
    eL = np.array([p.beta_Ls_err * np.cos(p.theta) for p in v])
    eG = np.array([p.beta_Gs_err * np.cos(p.theta) for p in v])
    mL, bL, smL, sbL, r2L = ols(x, yL)
    mG, smG = ols_origin(x, yG)
    ssr = np.sum((yG - mG * x) ** 2)
    sst = np.sum(yG**2)
    r.success = True
    r.slope, r.slope_err = mG, smG
    r.intercept, r.intercept_err = bL, sbL
    r.r2 = None
    r.extra.update(lorentz_r2=r2L, lorentz_slope=mL, lorentz_slope_err=smL, gauss_r2_origin=float(1 - ssr / sst) if sst > 0 else np.nan)
    _size_from_intercept(r, bL, sbL, K * lam)
    _strain_from_slope(r, mG, smG)
    r.plot = dict(x=list(x), y=list(yL * 1e3), yerr=list(eL * 1e3), y2=list(yG * 1e3), y2err=list(eG * 1e3),
                  xlabel="4 sin θ", ylabel="β cos θ (mrad)", scale=1e3,
                  fit1=(mL, bL), fit2=(mG, 0.0))
    return r


def halder_wagner(peaks, lam, K) -> MethodResult:
    r = MethodResult("hw", "Halder-Wagner", notes=NOTES["hw"])
    v = _valid(peaks)
    r.n = len(v)
    if len(v) < C.MIN_PEAKS_LINEAR:
        r.message = f"Needs ≥ {C.MIN_PEAKS_LINEAR} resolved reflections (have {len(v)})."
        return r
    bstar = np.array([p.beta_s * np.cos(p.theta) / lam for p in v])       # nm^-1
    bstar_err = np.array([p.beta_s_err * np.cos(p.theta) / lam for p in v])
    dstar = np.array([2 * np.sin(p.theta) / lam for p in v])               # nm^-1
    x = bstar / dstar**2                                                    # nm
    y = (bstar / dstar) ** 2                                                # dimensionless
    yerr = 2 * bstar * bstar_err / dstar**2
    m, c, sm, sc, r2 = ols(x, y)
    r.success = True
    r.slope, r.intercept, r.slope_err, r.intercept_err, r.r2 = m, c, sm, sc, r2
    # size: slope = K / D  (same significance rule as the W-H intercept)
    if m > 0 and (not np.isfinite(sm) or m > 2 * sm):
        r.size_nm = float(K / m)
        r.size_err = float(r.size_nm * sm / m) if np.isfinite(sm) else None
    else:
        top = m + 2 * sm if np.isfinite(sm) else m
        r.size_bound_nm = float(K / top) if top > 0 else None
        r.warnings.append("Slope not significantly above zero: size not resolved "
                          "(broadening is strain-dominated or the angular range is too short).")
    # strain: intercept = (2e)^2
    if c > 0 and (not np.isfinite(sc) or c > 2 * sc):
        e = np.sqrt(c) / 2
        r.strain, r.strain_err, r.strain_significant = float(e), float(sc / (8 * e)) if np.isfinite(sc) else None, True
    else:
        up = np.sqrt(max(c + 2 * (sc if np.isfinite(sc) else 0), 0)) / 2
        r.strain, r.strain_upper, r.strain_significant = max(0.0, float(np.sqrt(max(c, 0)) / 2)), float(up), False
    r.plot = dict(x=list(x), y=list(y * 1e3), yerr=list(yerr * 1e3),
                  xlabel="β*/d*² (nm)", ylabel="(β*/d*)² (×10⁻³)", scale=1e3)
    return r


def _size_from_intercept(r, b, sb, Klam):
    if b > 0 and (not np.isfinite(sb) or b > 2 * sb):
        r.size_nm = float(Klam / b)
        r.size_err = float(r.size_nm * sb / b) if np.isfinite(sb) else None
    else:
        top = b + 2 * sb if np.isfinite(sb) else b
        r.size_bound_nm = float(Klam / top) if top > 0 else None
        r.warnings.append("Intercept not significantly above zero: size not resolved "
                          "(broadening is strain-dominated or the angular range is too short).")


def _strain_from_slope(r, m, sm):
    if np.isfinite(sm) and abs(m) <= 2 * sm:
        r.strain_significant = False
        r.strain = float(m)
        r.strain_err = float(sm)
        r.strain_upper = float(max(m, 0) + 2 * sm)
    else:
        r.strain_significant = m > 0
        r.strain = float(m)
        r.strain_err = float(sm) if np.isfinite(sm) else None
        if m < 0:
            r.warnings.append("Negative slope: breadth falls with angle. This is not physical strain; "
                              "suspect anisotropic broadening, overlaps or a poor instrument profile.")
            r.strain_significant = False
            r.strain_upper = 0.0


METHODS = [
    ("scherrer", scherrer),
    ("single_line", single_line),
    ("monshi", monshi),
    ("wh", williamson_hall),
    ("voigt_wh", voigt_wh),
    ("hw", halder_wagner),
]


def run_methods(peaks: List[CorrectedPeak], lam: float, K: float) -> Dict[str, MethodResult]:
    out = {}
    for key, fn in METHODS:
        try:
            out[key] = fn(peaks, lam, K)
        except Exception as e:  # pragma: no cover - defensive
            out[key] = MethodResult(key, key, success=False, message=f"error: {e}")
    return out
