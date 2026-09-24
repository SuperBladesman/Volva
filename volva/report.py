"""
Reports and exports with full provenance.

Every export bundle contains:
  <name>_report.txt    human-readable report
  <name>_peaks.csv     per-reflection data
  <name>_summary.csv   one row per method
  <name>_session.json  everything needed to reproduce the result
"""

from __future__ import annotations

import csv
import json
import os
import platform

import numpy as np

from . import __version__
from .methods import _bound
from .pipeline import AnalysisResult

LEVEL_MARK = {"critical": "✗✗", "warning": "⚠", "info": "·", "ok": "✓"}


def _f(v, fmt="{:.4f}"):
    if v is None:
        return ""
    try:
        if not np.isfinite(v):
            return ""
    except TypeError:
        return str(v)
    return fmt.format(v)


def provenance(res: AnalysisResult) -> dict:
    import lmfit
    import scipy
    prof = res.profile
    return {
        "software": f"Völva {__version__}",
        "timestamp": res.timestamp,
        "python": platform.python_version(),
        "numpy": np.__version__, "scipy": scipy.__version__, "lmfit": lmfit.__version__,
        "platform": platform.platform(),
        "sample_file": res.pattern.path, "sample_sha256": res.pattern.sha256,
        "instrument_profile": None if prof is None else {
            "name": prof.name, "standard": prof.standard, "source_file": prof.source_file,
            "source_sha256": prof.source_sha256, "created": prof.created,
            "U": prof.U, "V": prof.V, "W": prof.W, "X": prof.X, "Y": prof.Y,
            "g_model": prof.g_model, "l_model": prof.l_model, "note": prof.note},
        "wavelength_nm": res.wavelength,
        "settings": res.settings.to_dict(),
        "K": res.settings.K,
        "background_half_width_deg": res.bg_half_width,
        "noise_sigma": res.noise,
    }


def _header_lines(res):
    pv = provenance(res)
    ip = pv["instrument_profile"]
    return [
        f"Völva {__version__} — crystallite size and strain analysis",
        f"Run: {res.timestamp}",
        f"Sample: {res.pattern.name}  (SHA-256 {res.pattern.sha256[:16]}…)",
        (f"Instrument profile: {ip['name']} from {ip['source_file']} (SHA-256 {ip['source_sha256'][:16]}…)"
         if ip else "Instrument profile: NONE — sizes are lower bounds"),
        f"λ(Kα1) = {res.wavelength['ka1']:.6f} nm; Kα2 {'fitted' if res.settings.doublet else 'not fitted'}; "
        f"K = {res.settings.K:.4g} ({res.settings.size_label}); breadth = integral breadth",
    ]


def text_report(res: AnalysisResult) -> str:
    L = []
    L += _header_lines(res)
    L.append("=" * 78)
    it = res.interpretation
    if it:
        L.append("RESULT")
        L.append("  " + it.headline)
        L.append("")
        L.append("WHAT THE METHODS SAY TOGETHER")
        for f in it.findings:
            L.append(f"  {LEVEL_MARK.get(f.level, '·')} {f.topic}: {f.text}")
    if res.warnings:
        L.append("")
        L.append("DATA / SETUP WARNINGS")
        for w in res.warnings:
            L.append(f"  ⚠ {w}")
    L.append("")
    L.append("METHOD COMPARISON")
    L.append(f"  {'Method':<30} {'Size':<26} {'Strain e':<40} {'R²':>6} {'n':>3}")
    for r in res.methods.values():
        r2 = f"{r.r2:.3f}" if r.r2 is not None and np.isfinite(r.r2) else ""
        if not r.success:
            L.append(f"  {r.name:<30} {'— ' + r.message}")
            continue
        note = "  (not used: profile shape undetermined)" if r.extra.get("shape_undetermined") else ""
        L.append(f"  {r.name:<30} {r.size_text():<26} {r.strain_text():<40} {r2:>6} {r.n:>3}{note}")
    L.append("  Size uncertainties: per-reflection methods give the SD between reflections; linear methods give")
    L.append("  the standard error of the regression. Upper limits are 2σ.")
    L.append("")
    L.append("REFLECTIONS")
    L.append(f"  {'2θ':>8} {'hkl':>6} {'use':>4} {'FWHM':>7} {'η':>5} {'G_s':>7} {'L_s':>7} "
             f"{'β/βi':>6} {'Scherrer':>10} {'flags'}")
    K, lam = res.settings.K, res.wavelength["ka1"]
    for p in res.peaks:
        D = K * lam / (p.beta_s * np.cos(p.theta)) if p.beta_s > 0 and not p.resolution_limited else None
        Bd = _bound(p, K, lam)
        Dt = f"{D:.1f}" if D else (f">{Bd:.0f}" if np.isfinite(Bd) else "—")
        L.append(f"  {p.two_theta:8.3f} {p.hkl:>6} {'yes' if p.use else 'no':>4} {p.fit.fwhm:7.4f} "
                 f"{p.fit.eta:5.2f} {p.G_s:7.4f} {p.L_s:7.4f} {p.ratio:6.2f} {Dt:>10} "
                 f"{', '.join(p.flags + p.fit.flags)}")
    L.append("")
    L.append("METHOD NOTES")
    for r in res.methods.values():
        L.append(f"  {r.name}: {r.notes}")
        for w in r.warnings:
            L.append(f"    ⚠ {w}")
    if res.profile is not None:
        L.append("")
        L.append(res.profile.summary())
    L.append("")
    L.append("LOG")
    L += ["  " + s for s in res.log]
    return "\n".join(L)


def write_peaks_csv(res: AnalysisResult, path: str):
    K, lam = res.settings.K, res.wavelength["ka1"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        for h in _header_lines(res):
            f.write(f"# {h}\n")
        w = csv.writer(f)
        w.writerow(["two_theta_deg", "two_theta_err_deg", "d_nm", "hkl", "use", "fwhm_deg", "fwhm_err_deg", "eta",
                    "eta_err", "skew", "integral_breadth_obs_deg", "G_obs_deg", "L_obs_deg", "G_inst_deg", "L_inst_deg",
                    "G_sample_deg", "L_sample_deg", "beta_sample_rad", "beta_sample_err_rad",
                    "beta_L_sample_rad", "beta_G_sample_rad", "broadening_ratio", "resolution_limited",
                    "scherrer_size_nm", "scherrer_lower_bound_nm", "group_R2", "group_size", "flags"])
        for p in res.peaks:
            rl = p.resolution_limited
            D = K * lam / (p.beta_s * np.cos(p.theta)) if p.beta_s > 0 and not rl else None
            B = _bound(p, K, lam) if rl else None
            w.writerow([_f(p.two_theta), _f(p.fit.center_err, "{:.5f}"), _f(p.d_nm, "{:.6f}"), p.hkl, p.use,
                        _f(p.fit.fwhm, "{:.5f}"), _f(p.fit.fwhm_err, "{:.5f}"), _f(p.fit.eta, "{:.3f}"),
                        _f(p.fit.eta_err, "{:.3f}"), _f(getattr(p.fit, "skew", 0.0), "{:.3f}"),
                        _f(p.beta_obs_deg, "{:.5f}"), _f(p.G_obs, "{:.5f}"),
                        _f(p.L_obs, "{:.5f}"), _f(p.G_inst, "{:.5f}"), _f(p.L_inst, "{:.5f}"),
                        _f(p.G_s, "{:.5f}"), _f(p.L_s, "{:.5f}"), _f(p.beta_s, "{:.6e}"),
                        _f(p.beta_s_err, "{:.3e}"), _f(p.beta_Ls, "{:.6e}"), _f(p.beta_Gs, "{:.6e}"),
                        _f(p.ratio, "{:.3f}"), rl, _f(D, "{:.2f}"), _f(B, "{:.0f}"),
                        _f(p.fit.r2), p.fit.group_size, ";".join(p.flags + p.fit.flags)])


def write_summary_csv(res: AnalysisResult, path: str):
    with open(path, "w", newline="", encoding="utf-8") as f:
        for h in _header_lines(res):
            f.write(f"# {h}\n")
        if res.interpretation:
            f.write(f"# Headline: {res.interpretation.headline}\n")
        w = csv.writer(f)
        w.writerow(["method", "success", "size_nm", "size_err_nm", "size_lower_bound_nm", "strain_percent",
                    "strain_err_percent", "strain_upper_limit_percent", "strain_resolved", "R2", "n_peaks",
                    "slope", "slope_err", "intercept", "intercept_err", "warnings", "message"])
        for r in res.methods.values():
            pc = lambda v: None if v is None else 100 * v  # noqa: E731
            w.writerow([r.name, r.success, _f(r.size_nm, "{:.2f}"), _f(r.size_err, "{:.2f}"),
                        _f(r.size_bound_nm, "{:.0f}"), _f(pc(r.strain), "{:.4f}"), _f(pc(r.strain_err), "{:.4f}"),
                        _f(pc(r.strain_upper), "{:.4f}"),
                        "" if r.strain_significant is None else r.strain_significant,
                        _f(r.r2), r.n, _f(r.slope, "{:.6e}"), _f(r.slope_err, "{:.3e}"),
                        _f(r.intercept, "{:.6e}"), _f(r.intercept_err, "{:.3e}"),
                        " | ".join(r.warnings), r.message])


def session_dict(res: AnalysisResult) -> dict:
    def clean(v):
        if isinstance(v, (np.floating, float)):
            return None if not np.isfinite(v) else float(v)
        if isinstance(v, (np.integer,)):
            return int(v)
        if isinstance(v, (np.bool_,)):
            return bool(v)
        if isinstance(v, np.ndarray):
            return [clean(x) for x in v.tolist()]
        if isinstance(v, dict):
            return {k: clean(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)):
            return [clean(x) for x in v]
        return v

    it = res.interpretation
    d = {
        "provenance": provenance(res),
        "headline": it.headline if it else "",
        "findings": [dict(level=f.level, topic=f.topic, text=f.text) for f in (it.findings if it else [])],
        "warnings": res.warnings,
        "methods": {k: {kk: vv for kk, vv in r.__dict__.items() if kk not in ("plot",)}
                    for k, r in res.methods.items()},
        "peaks": [dict(two_theta=p.two_theta, hkl=p.hkl, use=p.use, fwhm=p.fit.fwhm, fwhm_err=p.fit.fwhm_err,
                       eta=p.fit.eta, eta_err=p.fit.eta_err, skew=getattr(p.fit, "skew", 0.0), r2=p.fit.r2, G_obs=p.G_obs, L_obs=p.L_obs,
                       G_inst=p.G_inst, L_inst=p.L_inst, G_s=p.G_s, L_s=p.L_s, beta_s=p.beta_s,
                       beta_s_err=p.beta_s_err, beta_Ls=p.beta_Ls, beta_Gs=p.beta_Gs, ratio=p.ratio,
                       resolution_limited=p.resolution_limited, flags=p.flags + p.fit.flags)
                  for p in res.peaks],
        "rejected_fits": [dict(two_theta=f.center, reason=f.message, flags=f.flags)
                          for f in res.fits if not f.accepted],
        "instrument_profile_full": res.profile.to_dict() if res.profile else None,
        "log": res.log,
    }
    return clean(d)


def write_bundle(res: AnalysisResult, out_dir: str, base: str | None = None, figures: bool = True,
                 dpi: int = 300) -> list:
    os.makedirs(out_dir, exist_ok=True)
    base = base or os.path.splitext(res.pattern.name)[0]
    paths = []
    p = os.path.join(out_dir, f"{base}_report.txt")
    with open(p, "w", encoding="utf-8") as f:
        f.write(text_report(res))
    paths.append(p)
    p = os.path.join(out_dir, f"{base}_peaks.csv")
    write_peaks_csv(res, p)
    paths.append(p)
    p = os.path.join(out_dir, f"{base}_summary.csv")
    write_summary_csv(res, p)
    paths.append(p)
    p = os.path.join(out_dir, f"{base}_session.json")
    with open(p, "w", encoding="utf-8") as f:
        json.dump(session_dict(res), f, indent=2, ensure_ascii=False)
    paths.append(p)
    if figures:
        from matplotlib.figure import Figure
        from . import plots
        meta = {"Software": f"Völva {__version__}", "Source": res.pattern.name}
        for name, fn in [("pattern", lambda fg: plots.plot_pattern(fg, res)),
                         ("sizes", lambda fg: plots.plot_sizes(fg, res)),
                         ("comparison", lambda fg: plots.plot_comparison(fg, res)),
                         ("wh", lambda fg: plots.plot_linear(fg, res, "wh")),
                         ("voigt_wh", lambda fg: plots.plot_linear(fg, res, "voigt_wh")),
                         ("hw", lambda fg: plots.plot_linear(fg, res, "hw")),
                         ("instrument", lambda fg: plots.plot_instrument(fg, res.profile, res))]:
            fig = Figure(figsize=(9, 5.5))
            fn(fig)
            p = os.path.join(out_dir, f"{base}_{name}.png")
            fig.savefig(p, dpi=dpi, facecolor="white", metadata=meta)
            paths.append(p)
    return paths
