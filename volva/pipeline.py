"""
End-to-end analysis: pattern + instrument profile + settings -> results.
"""

from __future__ import annotations

import datetime as _dt
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional

import numpy as np

from . import __version__
from . import constants as C
from .background import auto_half_width, local_noise, noise_sigma, snip
from .correction import CorrectedPeak, correct_peaks
from .diagnostics import Interpretation, interpret
from .instrument import InstrumentProfile
from .io import Pattern
from .methods import MethodResult, run_methods
from .peaks import Candidate, PeakFit, detect_peaks, fit_all


@dataclass
class AnalysisSettings:
    wavelength: str = "Cu Ka"          # key of WAVELENGTHS, "From file" or "Custom"
    custom_ka1: float = 0.1540562
    custom_ka2: float = 0.1544390
    custom_ratio: float = 0.5
    doublet: bool = True
    size_convention: str = C.DEFAULT_SIZE_CONVENTION
    k_custom: Optional[float] = None
    bg_half_width: Optional[float] = None      # degrees; None = automatic
    prominence_sigma: float = C.DEFAULT_PROMINENCE_SIGMA
    min_rel_height: float = C.DEFAULT_MIN_REL_HEIGHT
    min_r2: float = C.DEFAULT_MIN_R2
    resolution_ratio: float = C.DEFAULT_RESOLUTION_RATIO
    two_theta_min: Optional[float] = None
    two_theta_max: Optional[float] = None
    mc_samples: int = C.DEFAULT_MC_SAMPLES
    seed: int = C.DEFAULT_SEED

    @property
    def K(self) -> float:
        if self.k_custom is not None:
            return float(self.k_custom)
        return C.SIZE_CONVENTIONS.get(self.size_convention, 1.0)

    @property
    def size_label(self) -> str:
        if self.k_custom is not None:
            return f"D (K = {self.k_custom:g})"
        return "D_sphere" if "Sphere" in self.size_convention else "L_vol"

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


def resolve_wavelength(settings: AnalysisSettings, pattern: Optional[Pattern] = None) -> dict:
    if settings.wavelength == "From file":
        if pattern is not None and pattern.wavelength:
            return dict(pattern.wavelength)
        return dict(C.WAVELENGTHS["Cu Ka"])
    if settings.wavelength == "Custom":
        return {"ka1": settings.custom_ka1, "ka2": settings.custom_ka2, "ratio": settings.custom_ratio}
    return dict(C.WAVELENGTHS[settings.wavelength])


@dataclass
class AnalysisResult:
    pattern: Pattern
    settings: AnalysisSettings
    wavelength: dict
    profile: Optional[InstrumentProfile]
    two_theta: np.ndarray = None
    intensity: np.ndarray = None
    background: np.ndarray = None
    bg_half_width: float = 0.0
    noise: float = 0.0
    candidates: List[Candidate] = field(default_factory=list)
    fits: List[PeakFit] = field(default_factory=list)
    peaks: List[CorrectedPeak] = field(default_factory=list)
    methods: Dict[str, MethodResult] = field(default_factory=dict)
    interpretation: Optional[Interpretation] = None
    log: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    critical: List[str] = field(default_factory=list)
    systematic: dict = field(default_factory=dict)
    timestamp: str = ""
    version: str = __version__


def _effective_doublet(doublet: bool, wl: dict) -> bool:
    return bool(doublet and wl.get("ratio", 0) > 0 and abs(wl.get("ka2", 0) - wl.get("ka1", 0)) > 1e-7)


def _compatibility(profile: Optional[InstrumentProfile], wl, settings, pattern, warnings, critical):
    if profile is None:
        return
    pw = profile.wavelength or {}
    if pw and abs(pw.get("ka1", 0) - wl["ka1"]) / wl["ka1"] > 2e-4:
        critical.append(
            f"The instrument profile was calibrated at λ = {pw.get('ka1'):.5f} nm but the sample uses "
            f"{wl['ka1']:.5f} nm. The profile does not apply to this sample.")
    d_prof = _effective_doublet(profile.doublet, pw or wl)
    d_samp = _effective_doublet(settings.doublet, wl)
    if d_prof != d_samp:
        critical.append(
            f"Kα2 handling differs: the standard was fitted {'with' if d_prof else 'without'} a Kα2 component, "
            f"the sample {'with' if d_samp else 'without'}. Widths are not comparable and a spurious strain or "
            "size error results. Use the same radiation setting for both, or recalibrate.")
    for w in profile.warnings:
        warnings.append(f"Instrument profile: {w}")


def run_analysis(pattern: Pattern, profile: Optional[InstrumentProfile],
                 settings: Optional[AnalysisSettings] = None, progress=None) -> AnalysisResult:
    settings = settings or AnalysisSettings()
    wl = resolve_wavelength(settings, pattern)
    res = AnalysisResult(pattern=pattern, settings=settings, wavelength=wl, profile=profile,
                         timestamp=_dt.datetime.now().isoformat(timespec="seconds"))
    log = res.log
    log.append(f"Völva {__version__} — {res.timestamp}")
    log.append(f"Sample: {pattern.name}  ({len(pattern.two_theta)} points, step {pattern.step:.4f}°)")
    for k, v in pattern.meta.items():
        if k.startswith("warning"):
            res.warnings.append(v)
    log.append(f"Wavelength: Kα1 {wl['ka1']:.6f} nm, Kα2 {wl['ka2']:.6f} nm, ratio {wl['ratio']:.2f}; "
               f"doublet {'fitted' if settings.doublet else 'not fitted'}")
    if settings.wavelength == "From file" and not pattern.wavelength:
        res.warnings.append("'From file' was selected but this file records no wavelength; Cu Kα was used.")
    _compatibility(profile, wl, settings, pattern, res.warnings, res.critical)

    tt, I = pattern.two_theta, pattern.intensity
    m = np.ones_like(tt, bool)
    if settings.two_theta_min is not None:
        m &= tt >= settings.two_theta_min
    if settings.two_theta_max is not None:
        m &= tt <= settings.two_theta_max
    tt, I = tt[m], I[m]
    res.two_theta, res.intensity = tt, I

    if progress:
        progress("Background")
    if settings.bg_half_width is None:
        hw, note = auto_half_width(tt, I)
        log.append(f"Background (SNIP): {note}")
    else:
        hw = float(settings.bg_half_width)
        log.append(f"Background (SNIP): half-width {hw:.2f}° (manual)")
    res.bg_half_width = hw
    res.background = snip(tt, I, hw)
    ysub = I - res.background
    sig_local = local_noise(I, res.background)
    res.noise = float(np.median(sig_local))
    step = float(np.median(np.diff(tt)))
    if profile is not None:
        # a sample peak cannot be much narrower than the instrument makes it
        min_fwhm = lambda t: max(C.MIN_FWHM_FRACTION_OF_INSTRUMENT * float(profile.total_fwhm_at(t)), 1.5 * step)
    else:
        min_fwhm = lambda t: 1.5 * step

    if progress:
        progress("Detecting peaks")
    res.candidates = detect_peaks(tt, ysub, wl, doublet=settings.doublet,
                                  prominence_sigma=settings.prominence_sigma,
                                  min_rel_height=settings.min_rel_height, noise=sig_local, log=log,
                                  min_fwhm=min_fwhm)
    log.append(f"Detected {len(res.candidates)} peaks (net area ≥ {settings.prominence_sigma:g}σ of the local "
               f"counting noise; median σ per point = {res.noise:.1f})")
    if progress:
        progress("Fitting profiles")
    res.fits = fit_all(tt, ysub, res.candidates, wl, doublet=settings.doublet,
                       min_r2=settings.min_r2, noise=sig_local, log=log, min_fwhm=min_fwhm)
    for f in res.fits:
        tag = "" if f.accepted else f"  REJECTED: {f.message}"
        fl = f" [{', '.join(f.flags)}]" if f.flags else ""
        log.append(f"  {f.center:8.3f}°  FWHM {f.fwhm:.4f}°  η {f.eta:.2f}  R² {f.r2:.4f}{fl}{tag}")
    widest = max((f.fwhm for f in res.fits if f.accepted), default=0.0)
    if widest > 0 and hw < 1.5 * widest:
        res.warnings.append(
            f"Background half-width ({hw:.2f}°) is less than 1.5 × the widest peak FWHM ({widest:.2f}°); "
            "the background may cut into peak tails. Increase it in Settings and compare.")

    if progress:
        progress("Instrument correction")
    res.peaks = correct_peaks(res.fits, wl, profile, resolution_ratio=settings.resolution_ratio,
                              n_mc=settings.mc_samples, seed=settings.seed)
    out = [p for p in res.peaks if "outside_calibrated_range" in p.flags]
    if out and profile is not None:
        lo, hi = profile.two_theta_range
        res.warnings.append(
            f"{len(out)} peak(s) lie outside the calibrated range {lo:.0f}–{hi:.0f}° "
            f"({', '.join(f'{p.two_theta:.1f}°' for p in out)}); the instrument profile is extrapolated there.")
    if progress:
        progress("Size/strain methods")
    _run(res)
    return res


def _run(res: AnalysisResult):
    s = res.settings
    res.methods = run_methods(res.peaks, res.wavelength["ka1"], s.K)
    _systematic(res)
    res.interpretation = interpret(res.methods, res.peaks, res.profile is not None, s.size_label,
                                   setup_critical=res.critical)
    it = res.interpretation
    r = res.methods.get(it.headline_key) if it.headline_key else None
    if r is not None and "size_sys" in r.extra and np.isfinite(it.headline_size):
        it.headline_size_sys = r.extra["size_sys"]
        extra = f" Instrument-calibration systematic: ± {r.extra['size_sys']:.1f} nm"
        if it.strain_state == "resolved" and "strain_sys" in r.extra:
            extra += f", ± {100 * r.extra['strain_sys']:.2g} % strain"
        it.headline += extra + "."


def _systematic(res: AnalysisResult):
    """
    Calibration uncertainty is common to all reflections, so it is not captured by
    the spread between reflections or by regression errors. Re-run the correction
    and every method with the instrument profile shifted by ±1σ; half the change is
    the systematic uncertainty. Because an error in the instrument's Gaussian width
    grows with tan θ it can mimic strain, so strain is re-tested against the
    combined statistical and systematic uncertainty.
    """
    res.systematic = {}
    if res.profile is None or not res.fits or not res.peaks:
        return
    from .correction import _ShiftedProfile
    use = {p.two_theta: (p.use, p.hkl) for p in res.peaks}
    runs = {}
    for mode in ("offset", "tilt"):
        for k in (-1.0, 1.0):
            peaks = correct_peaks(res.fits, res.wavelength, _ShiftedProfile(res.profile, k, mode),
                                  resolution_ratio=res.settings.resolution_ratio, n_mc=200, seed=res.settings.seed)
            for q in peaks:
                q.use, q.hkl = use.get(q.two_theta, (True, ""))
            runs[(mode, k)] = run_methods(peaks, res.wavelength["ka1"], res.settings.K)

    def half_diff(key, attr, mode):
        a, b = runs[(mode, -1.0)].get(key), runs[(mode, 1.0)].get(key)
        va, vb = (getattr(a, attr, None) if a else None), (getattr(b, attr, None) if b else None)
        if va is None or vb is None or not np.isfinite(va) or not np.isfinite(vb):
            return None
        return abs(vb - va) / 2

    for key, r in res.methods.items():
        s_sys = [v for v in (half_diff(key, "size_nm", m) for m in ("offset", "tilt")) if v is not None]
        e_sys = [v for v in (half_diff(key, "strain", m) for m in ("offset", "tilt")) if v is not None]
        if s_sys:
            r.extra["size_sys"] = float(max(s_sys))
        if e_sys:
            ss = float(max(e_sys))
            r.extra["strain_sys"] = ss
            if r.strain_significant and r.strain is not None:
                tot = np.sqrt((r.strain_err or 0.0) ** 2 + ss**2)
                if r.strain < 2 * tot:
                    r.strain_significant = False
                    r.strain_upper = float(max(r.strain, 0.0) + 2 * tot)
                    r.warnings.append(
                        f"Strain {100 * r.strain:.3f} % is within the instrument-calibration uncertainty "
                        f"(± {100 * ss:.3f} %); reported as an upper limit.")
        res.systematic[key] = {k2: r.extra[k2] for k2 in ("size_sys", "strain_sys") if k2 in r.extra}


def rerun_methods(res: AnalysisResult, use: Optional[List[bool]] = None,
                  hkl: Optional[List[str]] = None, settings: Optional[AnalysisSettings] = None):
    """Re-run the methods after the user changes peak selection, hkl labels or K."""
    if use is not None:
        for p, u in zip(res.peaks, use):
            p.use = bool(u)
    if hkl is not None:
        for p, h in zip(res.peaks, hkl):
            p.hkl = h
    if settings is not None:
        res.settings = settings
    _run(res)
    return res
