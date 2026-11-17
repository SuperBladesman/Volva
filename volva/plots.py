"""
Figures (matplotlib, no GUI dependency). Used by the GUI and the command line.

Colours: Okabe-Ito subset (colour-blind safe), always paired with a second
encoding (marker shape, line style or direct label).
"""

from __future__ import annotations

import numpy as np
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

from .methods import MethodResult

INK = "#1f1f1f"
MUTED = "#6b6b6b"
GRID = "#e6e6e6"
BLUE = "#0072B2"
VERM = "#D55E00"
GREEN = "#009E73"
ORANGE = "#E69F00"
PINK = "#CC79A7"
GREY = "#9a9a9a"


def _style(ax, xlabel=None, ylabel=None, title=None):
    ax.tick_params(colors=INK, labelsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#999999")
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    if xlabel:
        ax.set_xlabel(xlabel, color=INK, fontsize=10)
    if ylabel:
        ax.set_ylabel(ylabel, color=INK, fontsize=10)
    if title:
        ax.set_title(title, color=INK, fontsize=11, loc="left")


def _note(ax, text, loc="upper left"):
    x, ha = (0.02, "left") if "left" in loc else (0.98, "right")
    y, va = (0.97, "top") if "upper" in loc else (0.03, "bottom")
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=8, color=INK,
            bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#cccccc", alpha=0.9))


def _empty(fig, msg):
    fig.clear()
    ax = fig.add_subplot(111)
    ax.text(0.5, 0.5, msg, ha="center", va="center", color=MUTED, fontsize=12, transform=ax.transAxes)
    ax.set_axis_off()
    return fig


def _tight(fig):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fig.tight_layout()


# ---------------------------------------------------------------------------
def plot_raw(fig, pattern, background=None, title=None):
    fig.clear()
    ax = fig.add_subplot(111)
    ax.plot(pattern.two_theta, pattern.intensity, color=INK, lw=0.6, label="Data")
    if background is not None and len(background) == len(pattern.two_theta):
        ax.plot(pattern.two_theta, background, color=VERM, lw=1.2, ls="--", label="Background (SNIP)")
        ax.legend(frameon=False, fontsize=8)
    _style(ax, "2θ / °", "Intensity / counts", title or pattern.name)
    _tight(fig)
    return fig


def plot_pattern(fig, res):
    """Whole pattern: data, background, fitted groups, peak ticks and residual."""
    if res is None:
        return _empty(fig, "No analysis yet")
    fig.clear()
    gs = GridSpec(2, 1, height_ratios=[4, 1], hspace=0.05, figure=fig)
    ax = fig.add_subplot(gs[0])
    axr = fig.add_subplot(gs[1], sharex=ax)
    tt, I, bg = res.two_theta, res.intensity, res.background
    ax.plot(tt, I, color=INK, lw=0.5, label="Data")
    ax.plot(tt, bg, color=GREY, lw=1.0, ls="--", label=f"Background (SNIP ±{res.bg_half_width:.2f}°)")
    seen = set()
    for f in res.fits:
        if f.fit_x is None or f.group_id in seen:
            continue
        seen.add(f.group_id)
        ax.plot(f.fit_x, f.fit_total + np.interp(f.fit_x, tt, bg), color=VERM, lw=1.2)
    use = {id(p.fit): p for p in res.peaks}
    ytop = I.max()
    for f in res.fits:
        p = use.get(id(f))
        if not f.accepted:
            ax.plot(f.center, ytop * 1.04, marker="x", color=GREY, ms=6, mew=1.2)
        elif p is not None and p.resolution_limited:
            ax.plot(f.center, ytop * 1.04, marker="^", color=ORANGE, ms=6, mec=INK, mew=0.4)
        elif p is not None and not p.use:
            ax.plot(f.center, ytop * 1.04, marker="|", color=GREY, ms=10, mew=1.5)
        else:
            ax.plot(f.center, ytop * 1.04, marker="|", color=BLUE, ms=10, mew=1.5)
    handles = [Line2D([], [], color=INK, lw=0.8, label="Data"),
               Line2D([], [], color=GREY, lw=1, ls="--", label="Background"),
               Line2D([], [], color=VERM, lw=1.2, label="Fitted profile"),
               Line2D([], [], color=BLUE, marker="|", ls="", ms=10, mew=1.5, label="Peak used"),
               Line2D([], [], color=ORANGE, marker="^", ls="", ms=6, label="Resolution-limited"),
               Line2D([], [], color=GREY, marker="x", ls="", ms=6, label="Rejected / unused")]
    ax.legend(handles=handles, frameon=False, fontsize=8, ncol=3, loc="upper right")
    ax.set_ylim(top=ytop * 1.25)
    _style(ax, None, "Intensity / counts", res.pattern.name)
    ax.tick_params(labelbottom=False)
    for f in res.fits:
        if f.residual is not None:
            axr.plot(f.data_x, f.residual, color=MUTED, lw=0.5)
    axr.axhline(0, color=INK, lw=0.5)
    _style(axr, "2θ / °", "Resid.")
    _tight(fig)
    return fig


def plot_peak(fig, res, index):
    """One peak: data, group fit, Kα1/Kα2 components, residual and numbers."""
    if res is None or not res.peaks:
        return _empty(fig, "No fitted peaks")
    index = max(0, min(index, len(res.peaks) - 1))
    p = res.peaks[index]
    f = p.fit
    fig.clear()
    gs = GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05, figure=fig)
    ax = fig.add_subplot(gs[0])
    axr = fig.add_subplot(gs[1], sharex=ax)
    ax.plot(f.data_x, f.data_y, "o", color=BLUE, ms=3, alpha=0.7, label="Data − background")
    ax.plot(f.fit_x, f.fit_total, color=VERM, lw=1.6, label="Group fit")
    ax.plot(f.fit_x, f.fit_peak_ka1 + f.fit_bg, color=GREEN, lw=1.0, ls="--", label="Kα1")
    if f.doublet:
        ax.plot(f.fit_x, f.fit_peak - f.fit_peak_ka1 + f.fit_bg, color=ORANGE, lw=1.0, ls=":", label="Kα2")
    ax.plot(f.fit_x, f.fit_bg, color=GREY, lw=0.8, ls="-.", label="Local background")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    K, lam = res.settings.K, res.wavelength["ka1"]
    D = K * lam / (p.beta_s * np.cos(p.theta)) if p.beta_s > 0 else np.nan
    lines = [
        f"2θ = {p.two_theta:.3f}°   d = {10 * p.d_nm:.4f} Å" + (f"   hkl {p.hkl}" if p.hkl else ""),
        f"FWHM = {f.fwhm:.4f}° ± {f.fwhm_err:.4f}   η = {f.eta:.2f} ± {f.eta_err:.2f}"
        + (f"   skew {f.skew:+.2f}" if getattr(f, "skew", 0.0) else ""),
        f"Observed  G {p.G_obs:.4f}°  L {p.L_obs:.4f}°",
        f"Instrument G {p.G_inst:.4f}°  L {p.L_inst:.4f}°",
        f"Sample    G {p.G_s:.4f}°  L {p.L_s:.4f}°",
        f"β_obs/β_inst = {p.ratio:.2f}   group R² = {f.r2:.4f}",
    ]
    if p.resolution_limited:
        lines.append("RESOLUTION-LIMITED: size not measurable")
    elif np.isfinite(D):
        lines.append(f"Scherrer {res.settings.size_label} = {D:.1f} nm")
    if f.flags:
        lines.append("flags: " + ", ".join(f.flags))
    _note(ax, "\n".join(lines))
    _style(ax, None, "Intensity", f"Peak {index + 1} of {len(res.peaks)}")
    ax.tick_params(labelbottom=False)
    axr.plot(f.data_x, f.residual, color=MUTED, lw=0.8)
    axr.axhline(0, color=INK, lw=0.5)
    if res.noise:
        for s in (-3, 3):
            axr.axhline(s * res.noise, color=VERM, lw=0.5, ls="--")
    _style(axr, "2θ / °", "Resid.")
    _tight(fig)
    return fig


def plot_gallery(fig, res):
    if res is None or not res.peaks:
        return _empty(fig, "No fitted peaks")
    fig.clear()
    n = len(res.peaks)
    nc = int(np.ceil(np.sqrt(n * 1.5)))
    nr = int(np.ceil(n / nc))
    for i, p in enumerate(res.peaks):
        ax = fig.add_subplot(nr, nc, i + 1)
        f = p.fit
        ax.plot(f.data_x, f.data_y, ".", color=BLUE, ms=2)
        ax.plot(f.fit_x, f.fit_total, color=VERM, lw=1)
        tag = " RL" if p.resolution_limited else ""
        tag += " ×" if not p.use else ""
        ax.set_title(f"{p.two_theta:.2f}°{tag}  R² {f.r2:.3f}", fontsize=7, color=INK)
        ax.tick_params(labelsize=6)
        ax.set_yticks([])
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
    _tight(fig)
    return fig


def plot_instrument(fig, profile, res=None):
    """
    Instrument profile: total FWHM (Caglioti) with the derived G and L components, the
    Lorentzian fraction η and peak asymmetry, and (after an analysis) sample / instrument breadth.
    """
    from .profiles import compose
    if profile is None:
        return _empty(fig, "No instrument profile loaded")
    fig.clear()
    ncol = 3 if (res is not None and res.peaks) else 2
    gs = GridSpec(1, ncol, width_ratios=[3, 2, 2][:ncol], wspace=0.32, figure=fig)
    lo, hi = profile.two_theta_range
    x = np.linspace(max(5, lo - 10), min(165, hi + 10), 300)
    pk = profile.peaks

    def arr(key, pts):
        return np.array([q.get(key, np.nan) for q in pts], float)

    def fwhm_eta(pts):
        if pts and "fwhm" in pts[0]:
            return arr("fwhm", pts), arr("eta", pts)
        if not pts:
            return np.array([]), np.array([])
        F, E = compose(arr("G", pts), arr("L", pts))
        return np.atleast_1d(F), np.atleast_1d(E)

    used = [q for q in pk if q.get("used", True)]
    out = [q for q in pk if not q.get("used", True)]
    Fu, Eu = fwhm_eta(used)
    Fo, Eo = fwhm_eta(out)

    # --- widths ----------------------------------------------------------------
    ax = fig.add_subplot(gs[0])
    ax.axvspan(lo, hi, color="#f3f3f3", zorder=0)
    total = (profile.fwhm_at(x) if profile.profile_model == "total"
             else np.atleast_1d(compose(profile.G_at(x), profile.L_at(x))[0]))
    ax.plot(x, total, color=INK, lw=1.8, label="Total FWHM (model)")
    ax.plot(x, profile.G_at(x), color=BLUE, lw=1.3, label="Gaussian part G")
    ax.plot(x, profile.L_at(x), color=VERM, lw=1.3, ls="--", label="Lorentzian part L")
    ax.plot(arr("two_theta", used), Fu, "o", color=INK, ms=5, mfc="white", label="Standard peaks")
    if out:
        ax.plot(arr("two_theta", out), Fo, "x", color=GREY, ms=6, label="Excluded")
    if res is not None and res.peaks:
        Fs = np.array([compose(p.G_obs, p.L_obs)[0] for p in res.peaks], float).ravel()
        ax.plot([p.two_theta for p in res.peaks], Fs, "^", color=GREEN, ms=5, label="Sample FWHM")
    ax.set_ylim(bottom=0)
    ax.legend(frameon=False, fontsize=7)
    _style(ax, "2θ / °", "FWHM / °", "Instrument width (shaded: calibrated range)")

    # --- shape ------------------------------------------------------------------
    ax2 = fig.add_subplot(gs[1])
    ax2.axvspan(lo, hi, color="#f3f3f3", zorder=0)
    ee = arr("eta_err", used) if used and "eta_err" in used[0] else None
    ax2.errorbar(arr("two_theta", used), Eu, yerr=ee, fmt="o", color=BLUE, ms=5, capsize=2, lw=1,
                 label="η (standard)")
    Ex = (profile.eta_at(x) if profile.profile_model == "total"
          else np.atleast_1d(compose(profile.G_at(x), profile.L_at(x))[1]))
    ax2.plot(x, Ex, color=BLUE, lw=1.5, label="η (model)")
    if used and "skew" in used[0]:
        ax2.plot(arr("two_theta", used), arr("skew", used), "s", color=ORANGE, ms=5, mec=INK, mew=0.4,
                 label="Asymmetry (skew)")
    ax2.axhline(0, color=INK, lw=0.5)
    ax2.set_ylim(-0.1, 1.1)
    ax2.legend(frameon=False, fontsize=7, loc="lower left")
    _style(ax2, "2θ / °", "η and skew", "Peak shape")

    # --- sample / instrument ------------------------------------------------------
    if ncol == 3:
        ax3 = fig.add_subplot(gs[2])
        thr = res.settings.resolution_ratio
        xs = [p.two_theta for p in res.peaks]
        ys = [p.ratio for p in res.peaks]
        cols = [ORANGE if (p.resolution_limited or "near_resolution_limit" in p.flags) else BLUE for p in res.peaks]
        ax3.scatter(xs, ys, c=cols, s=30, edgecolors=INK, linewidths=0.4)
        ax3.axhline(thr, color=VERM, ls="--", lw=1)
        ax3.text(0.98, thr, f" limit {thr:g}", transform=ax3.get_yaxis_transform(), ha="right", va="bottom",
                 fontsize=8, color=VERM)
        ax3.set_ylim(bottom=0.8)
        _style(ax3, "2θ / °", "β_obs/β_inst", "Sample / instrument breadth")
    _tight(fig)
    return fig


def plot_sizes(fig, res):
    """Per-reflection sizes (Scherrer, single-line) and single-line strains vs 2θ."""
    if res is None or not res.methods:
        return _empty(fig, "No results")
    sch = res.methods.get("scherrer")
    sl = res.methods.get("single_line")
    fig.clear()
    gs = GridSpec(2, 1, height_ratios=[3, 2], hspace=0.12, figure=fig)
    ax = fig.add_subplot(gs[0])
    lab = res.settings.size_label
    if sch and sch.plot:
        ax.errorbar(sch.plot["x"], sch.plot["y"], yerr=sch.plot["yerr"], fmt="o", color=BLUE, ms=6,
                    capsize=2, lw=1, label="Scherrer (all breadth = size)")
        if sch.plot.get("xb"):
            yb = [v for v in sch.plot["yb"]]
            ax.plot(sch.plot["xb"], yb, "^", color=ORANGE, ms=8, mec=INK, mew=0.4,
                    label="Resolution-limited (lower bound)")
        if sch.size_nm:
            ax.axhline(sch.size_nm, color=BLUE, lw=0.8, ls=":")
    if sl and sl.plot.get("x"):
        ax.errorbar(np.array(sl.plot["x"]) + 0.3, sl.plot["y"], yerr=sl.plot["yerr"], fmt="s", color=VERM,
                    ms=5, mfc="white", capsize=2, lw=1, label="Single-line (Lorentzian part)")
    for p in res.peaks:
        if p.hkl and p.use:
            ax.annotate(p.hkl, (p.two_theta, 0), xycoords=("data", "axes fraction"), xytext=(0, 3),
                        textcoords="offset points", ha="center", fontsize=7, color=MUTED)
    ax.legend(frameon=False, fontsize=8)
    ax.set_ylim(bottom=0)
    _style(ax, None, f"{lab} / nm", "Size by reflection")
    ax.tick_params(labelbottom=False)
    ax2 = fig.add_subplot(gs[1], sharex=ax)
    if sl and sl.plot.get("xs"):
        ax2.errorbar(sl.plot["xs"], sl.plot["ys"], yerr=sl.plot["yserr"], fmt="D", color=GREEN, ms=5,
                     capsize=2, lw=1, label="Single-line strain (Gaussian part)")
        ax2.axhline(0, color=INK, lw=0.5)
        ax2.legend(frameon=False, fontsize=8)
    _style(ax2, "2θ / °", "Strain e / %")
    _tight(fig)
    return fig


_LIN_TITLES = {"monshi": "Monshi-Scherrer", "wh": "Williamson-Hall (UDM)",
               "voigt_wh": "Voigt Williamson-Hall", "hw": "Halder-Wagner"}


def _band(ax, x, y, m, b, color):
    """Line with 95% confidence band from OLS."""
    from scipy import stats
    x = np.asarray(x)
    y = np.asarray(y)
    n = len(x)
    xs = np.linspace(min(0, x.min()), x.max() * 1.05, 200)
    ys = m * xs + b
    ax.plot(xs, ys, color=color, lw=1.4)
    if n > 2:
        res = y - (m * x + b)
        s = np.sqrt(np.sum(res**2) / (n - 2))
        sxx = np.sum((x - x.mean()) ** 2)
        t = stats.t.ppf(0.975, n - 2)
        se = s * np.sqrt(1 / n + (xs - x.mean()) ** 2 / sxx)
        ax.fill_between(xs, ys - t * se, ys + t * se, color=color, alpha=0.12, lw=0)


def plot_linear(fig, res, key):
    if res is None or key not in res.methods:
        return _empty(fig, "No results")
    r: MethodResult = res.methods[key]
    if not r.success or not r.plot:
        return _empty(fig, f"{_LIN_TITLES.get(key, key)}: {r.message or 'not available'}")
    fig.clear()
    ax = fig.add_subplot(111)
    x = np.array(r.plot["x"])
    y = np.array(r.plot["y"])
    sc = r.plot.get("scale", 1.0)
    if key == "voigt_wh":
        y2 = np.array(r.plot["y2"])
        ax.errorbar(x, y, yerr=r.plot["yerr"], fmt="o", color=VERM, ms=6, capsize=2, label="Lorentzian β_L cos θ (size)")
        ax.errorbar(x, y2, yerr=r.plot["y2err"], fmt="s", color=BLUE, ms=6, mfc="white", capsize=2,
                    label="Gaussian β_G cos θ (strain)")
        mL, bL = r.plot["fit1"]
        _band(ax, x, y, mL * sc, bL * sc, VERM)
        xs = np.linspace(0, x.max() * 1.05, 50)
        ax.plot(xs, r.plot["fit2"][0] * sc * xs, color=BLUE, lw=1.4, ls="--")
        ax.legend(frameon=False, fontsize=8, loc="lower right")
    else:
        ax.errorbar(x, y, yerr=r.plot.get("yerr"), fmt="o", color=BLUE, ms=6, capsize=2)
        _band(ax, x, y, r.slope * sc, r.intercept * sc, VERM)
        for p_x, p_y, p in zip(x, y, [q for q in res.peaks if q.use and not q.resolution_limited and q.beta_s > 0]):
            if p.hkl:
                ax.annotate(p.hkl, (p_x, p_y), xytext=(4, 4), textcoords="offset points", fontsize=7, color=MUTED)
    lines = [f"n = {r.n}"]
    if r.r2 is not None and np.isfinite(r.r2) and key != "voigt_wh":
        lines.append(f"R² = {r.r2:.3f}")
    lines.append(f"Size: {r.size_text()}")
    if key != "monshi":
        lines.append(f"Strain: {r.strain_text()}")
    else:
        lines.append(f"Slope = {r.slope:.3f} ± {r.slope_err:.3f} (1 for pure size)")
    for w in r.warnings[:2]:
        lines.append("⚠ " + w[:70] + ("…" if len(w) > 70 else ""))
    _note(ax, "\n".join(lines))
    _style(ax, r.plot["xlabel"], r.plot["ylabel"], _LIN_TITLES.get(key, r.name))
    _tight(fig)
    return fig


def plot_comparison(fig, res):
    """Dot-and-whisker comparison of size and strain across methods."""
    if res is None or not res.methods:
        return _empty(fig, "No results")
    order = ["scherrer", "single_line", "monshi", "wh", "voigt_wh", "hw"]
    rs = [res.methods[k] for k in order if k in res.methods]
    fig.clear()
    gs = GridSpec(1, 2, wspace=0.08, figure=fig)
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharey=ax)
    ypos = np.arange(len(rs))[::-1]
    head = res.interpretation.headline_size if res.interpretation else np.nan
    for yv, r in zip(ypos, rs):
        if r.size_nm is not None and np.isfinite(r.size_nm):
            ax.errorbar(r.size_nm, yv, xerr=r.size_err if r.size_err else None, fmt="o", color=BLUE, ms=7,
                        capsize=3, lw=1.2)
            ax.annotate(f"{r.size_nm:.1f}", (r.size_nm, yv), xytext=(0, 7), textcoords="offset points",
                        ha="center", fontsize=8, color=INK)
        elif r.size_bound_nm:
            ax.plot(r.size_bound_nm, yv, marker=">", color=ORANGE, ms=8, mec=INK, mew=0.4)
            ax.annotate(f"> {r.size_bound_nm:.0f}", (r.size_bound_nm, yv), xytext=(0, 7),
                        textcoords="offset points", ha="center", fontsize=8, color=INK)
        else:
            ax.annotate("not available", (0.02, yv), xycoords=("axes fraction", "data"), fontsize=8, color=MUTED,
                        va="center")
    if np.isfinite(head):
        ax.axvline(head, color=VERM, lw=0.8, ls="--")
    short = {"scherrer": "Scherrer", "single_line": "Single-line Voigt", "monshi": "Monshi",
             "wh": "W-H (UDM)", "voigt_wh": "Voigt W-H", "hw": "Halder-Wagner"}
    ax.set_yticks(ypos)
    ax.set_yticklabels([short.get(r.key, r.name) for r in rs], fontsize=9)
    it = res.interpretation
    if it is not None and it.strain_state == "resolved":
        for yv, r in zip(ypos, rs):
            if r.key == "scherrer" and r.size_nm:
                ax.annotate("lower bound (strain present)", (r.size_nm, yv), xytext=(0, -12),
                            textcoords="offset points", ha="center", fontsize=7, color=MUTED)
    ax.set_xlim(left=0)
    _style(ax, f"{res.settings.size_label} / nm", None, "Size")
    for yv, r in zip(ypos, rs):
        if r.strain_significant is False and r.strain_upper is not None:
            ax2.plot([0, 100 * r.strain_upper], [yv, yv], color=GREY, lw=3, solid_capstyle="round")
            ax2.annotate(f"< {100 * r.strain_upper:.3f}", (100 * r.strain_upper, yv), xytext=(4, 0),
                         textcoords="offset points", va="center", fontsize=8, color=MUTED)
        elif r.strain is not None and np.isfinite(r.strain):
            ax2.errorbar(100 * r.strain, yv, xerr=100 * r.strain_err if r.strain_err else None, fmt="s",
                         color=GREEN, ms=7, capsize=3, lw=1.2)
            ax2.annotate(f"{100 * r.strain:.3f}", (100 * r.strain, yv), xytext=(0, 7), textcoords="offset points",
                         ha="center", fontsize=8, color=INK)
    ax2.axvline(0, color=INK, lw=0.5)
    ax2.tick_params(labelleft=False)
    _style(ax2, "Strain e / %", None, "Strain (grey bar = upper limit)")
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.12, top=0.9, wspace=0.08)
    return fig
