"""
Peak detection and profile fitting.

* Detection on background-subtracted, lightly smoothed data with a
  noise-referenced prominence threshold.
* Resolved Kα2 components are recognised and not treated as separate peaks.
* Overlapping peaks are fitted together (grouped fitting).
* Each peak: pseudo-Voigt Kα1 + constrained Kα2 (exact Bragg-law position,
  fixed intensity ratio, same width and shape) on a local linear background.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
from lmfit.lineshapes import pvoigt
from scipy.signal import find_peaks, savgol_filter

from . import constants as C
from .background import noise_sigma
from .profiles import pv_integral_breadth


def split_pvoigt(x, amplitude=1.0, center=0.0, sigma=0.1, fraction=0.5, skew=0.0):
    """
    Asymmetric (split) pseudo-Voigt in lmfit's convention.

    FWHM = 2*sigma and area = amplitude, exactly as lmfit's PseudoVoigtModel.
    The low-angle half-width is sigma*(1+skew), the high-angle half-width
    sigma*(1-skew), so skew > 0 gives the low-angle tail typical of axial
    divergence. Height, FWHM and integral breadth do not depend on skew, so the
    breadth mathematics is unchanged.
    """
    s = np.where(x < center, 1.0 + skew, 1.0 - skew)
    return pvoigt((x - center) / s, amplitude, 0.0, sigma, fraction)


@dataclass
class Candidate:
    pos: float
    height: float
    fwhm_est: float
    prominence: float
    snr: float


@dataclass
class PeakFit:
    center: float
    fwhm: float
    eta: float
    amplitude: float
    height: float
    center_err: float = np.nan
    fwhm_err: float = np.nan
    eta_err: float = np.nan
    skew: float = 0.0
    snr_fit: float = np.nan        # fitted area / its standard error
    chi2_red: float = np.nan       # group reduced chi-square against the local counting noise
    cov_fwhm_eta: Optional[np.ndarray] = None
    r2: float = 0.0
    snr: float = 0.0
    group_id: int = -1
    group_size: int = 1
    doublet: bool = False
    accepted: bool = False
    flags: List[str] = field(default_factory=list)
    message: str = ""
    # plotting data
    data_x: Optional[np.ndarray] = None
    data_y: Optional[np.ndarray] = None
    fit_x: Optional[np.ndarray] = None
    fit_total: Optional[np.ndarray] = None
    fit_peak: Optional[np.ndarray] = None
    fit_peak_ka1: Optional[np.ndarray] = None
    fit_bg: Optional[np.ndarray] = None
    residual: Optional[np.ndarray] = None

    @property
    def integral_breadth(self) -> float:
        return float(pv_integral_breadth(self.fwhm, self.eta))


def ka2_position(two_theta_ka1, wl):
    """Exact Kα2 position (degrees) from Bragg's law."""
    r = wl["ka2"] / wl["ka1"]
    s = np.clip(r * np.sin(np.deg2rad(np.asarray(two_theta_ka1) / 2.0)), -1, 1)
    return 2.0 * np.rad2deg(np.arcsin(s))


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------
def detect_peaks(two_theta, y_sub, wl, doublet=True,
                 prominence_sigma=C.DEFAULT_PROMINENCE_SIGMA,
                 min_rel_height=C.DEFAULT_MIN_REL_HEIGHT,
                 noise=None, log=None, min_fwhm=None) -> List[Candidate]:
    step = float(np.median(np.diff(two_theta)))
    w = max(5, int(round(0.06 / step)) | 1)
    w = min(w, (len(y_sub) // 2) * 2 - 1)
    ys = savgol_filter(y_sub, w, 2) if w >= 5 else y_sub.copy()
    sig = noise if noise is not None else noise_sigma(y_sub)
    sig_arr = np.broadcast_to(np.asarray(sig, float), y_sub.shape)   # scalar or point-by-point
    ymax = float(np.max(ys))
    idx, props = find_peaks(
        ys,
        prominence=float(np.min(sig_arr)),
        height=min_rel_height * ymax,
        distance=max(2, int(round(0.03 / step))),
        width=1.5,
    )
    # Significance of each feature = its net counts over +/- 1 FWHM above a LOCAL straight
    # baseline (drawn between points 2 FWHM either side), divided by the counting noise of
    # that net area. Unlike a height or prominence test this rejects the many small noise
    # excursions of a long scan, including those sitting on the flanks of broad peaks or on
    # a curved background, and it uses the local noise, which is larger where the
    # background is higher.
    n = len(y_sub)
    cands = []
    n_weak = n_narrow = 0
    for k, i in enumerate(idx):
        fw = float(props["widths"][k] * step)
        mf = float(min_fwhm(two_theta[i])) if min_fwhm is not None else 2 * step
        if min_fwhm is not None and fw < 0.5 * mf:
            n_narrow += 1
            continue
        nh = max(2, int(np.ceil(max(fw, mf, 2 * step) / step)))
        a0, a1 = max(0, i - nh), min(n, i + nh + 1)
        e0, e1 = max(0, i - 2 * nh), min(n - 1, i + 2 * nh)
        s0, s1 = slice(max(0, e0 - 2), e0 + 3), slice(max(0, e1 - 2), min(n, e1 + 3))
        b0, b1 = float(np.mean(y_sub[s0])), float(np.mean(y_sub[s1]))
        xs = np.arange(a0, a1)
        base = b0 + (b1 - b0) * (xs - e0) / max(e1 - e0, 1)
        area = float(np.sum(y_sub[a0:a1] - base))
        nw = a1 - a0
        var = float(np.sum(sig_arr[a0:a1] ** 2)) + (nw / 2.0) ** 2 * (
            np.mean(sig_arr[s0] ** 2) / max(len(y_sub[s0]), 1) + np.mean(sig_arr[s1] ** 2) / max(len(y_sub[s1]), 1))
        snr = area / np.sqrt(var) if var > 0 else 0.0
        if snr < prominence_sigma:
            n_weak += 1
            continue
        cands.append(Candidate(pos=float(two_theta[i]), height=float(ys[i]),
                               fwhm_est=max(fw, 2 * step),
                               prominence=float(props["prominences"][k]), snr=snr))
    if log is not None and (n_weak or n_narrow):
        log.append(f"  Detection: {n_weak} feature(s) below the local noise threshold and {n_narrow} "
                   "narrower than the instrument allows were ignored")
    cands = _merge_split(cands, log)
    if doublet and wl.get("ratio", 0) > 0 and wl["ka2"] != wl["ka1"]:
        cands = _drop_ka2(cands, wl, step, log)
    return cands


def _merge_split(cands, log):
    """A noisy top can split one broad peak into two maxima; keep the more significant one."""
    cands = sorted(cands, key=lambda c: c.pos)
    out = []
    for c in cands:
        if out and abs(c.pos - out[-1].pos) < 0.5 * max(c.fwhm_est, out[-1].fwhm_est):
            weak, strong = (c, out[-1]) if c.snr < out[-1].snr else (out[-1], c)
            if log is not None:
                log.append(f"  {weak.pos:8.3f}°: part of the peak at {strong.pos:.3f}° (split maximum), merged")
            out[-1] = strong
        else:
            out.append(c)
    return out


def _drop_ka2(cands, wl, step, log):
    """Remove detections that are the resolved Kα2 partner of a stronger peak."""
    keep = [True] * len(cands)
    for i, ci in enumerate(cands):
        p2 = float(ka2_position(ci.pos, wl))
        split = p2 - ci.pos
        tol = max(0.35 * split, 2 * step)
        for j, cj in enumerate(cands):
            if j == i or not keep[j]:
                continue
            if abs(cj.pos - p2) < tol and cj.height < 0.85 * ci.height and cj.pos > ci.pos:
                keep[j] = False
                if log is not None:
                    log.append(f"  {cj.pos:8.3f}°: recognised as Kα2 of {ci.pos:.3f}° (not a separate peak)")
    return [c for c, k in zip(cands, keep) if k]


# ---------------------------------------------------------------------------
# Grouping
# ---------------------------------------------------------------------------
def group_candidates(cands: List[Candidate], wl, doublet=True, reach=2.5, max_group=8):
    """Group candidates whose fitting windows overlap."""
    if not cands:
        return []
    cands = sorted(cands, key=lambda c: c.pos)

    def span(c):
        d = float(ka2_position(c.pos, wl) - c.pos) if doublet else 0.0
        return c.pos - reach * c.fwhm_est, c.pos + reach * c.fwhm_est + d

    groups = [[cands[0]]]
    for c in cands[1:]:
        lo, _ = span(c)
        _, hi_prev = span(groups[-1][-1])
        if lo < hi_prev:
            groups[-1].append(c)
        else:
            groups.append([c])
    out = []
    for g in groups:
        while len(g) > max_group:
            gaps = np.diff([c.pos for c in g])
            k = int(np.argmax(gaps[1:-1])) + 1 if len(gaps) > 2 else len(g) // 2
            out.append(g[:k + 1])
            g = g[k + 1:]
        out.append(g)
    return out


# ---------------------------------------------------------------------------
# Fitting
#
# Each group is fitted with scipy.optimize.least_squares (trust-region
# reflective, native bounds) on a flat parameter vector. The Kα2 component is
# computed inside the model from the Kα1 parameters (exact Bragg position,
# fixed intensity ratio, same width, shape and asymmetry), so there are no
# constraint expressions to evaluate. This is 10-50x faster than the
# equivalent lmfit composite model and gives the same optimum.
# ---------------------------------------------------------------------------
_S2L2 = np.sqrt(2.0 * np.log(2.0))
_NP = 5                      # per-peak parameters: amplitude, center, sigma, fraction, skew


def _spv(x, amp, cen, sig, eta, skew):
    """split_pvoigt, inlined for speed (same definition)."""
    u = (x - cen) / np.where(x < cen, 1.0 + skew, 1.0 - skew)
    sg = sig / _S2L2
    g = np.exp(-0.5 * (u / sg) ** 2) / (sg * np.sqrt(2.0 * np.pi))
    lo = (sig / np.pi) / (u * u + sig * sig)
    return amp * ((1.0 - eta) * g + eta * lo)


def _ka2_center(cen, r):
    return 2.0 * np.degrees(np.arcsin(np.clip(r * np.sin(np.radians(cen) / 2.0), -1.0, 1.0)))


class _GroupModel:
    """Linear background + n split pseudo-Voigts (+ tied Kα2) on a flat vector."""

    def __init__(self, n, use_d, r, ratio):
        self.n, self.use_d, self.r, self.ratio = n, use_d, r, ratio

    def components(self, x, v):
        out = []
        for k in range(self.n):
            a, c, sg, e, sk = v[2 + _NP * k: 2 + _NP * (k + 1)]
            p1 = _spv(x, a, c, sg, e, sk)
            p2 = _spv(x, self.ratio * a, _ka2_center(c, self.r), sg, e, sk) if self.use_d else 0.0 * x
            out.append((p1, p2))
        return out

    def __call__(self, x, v):
        y = v[0] + v[1] * x
        for p1, p2 in self.components(x, v):
            y = y + p1 + p2
        return y


def _solve(model, x, y, v0, lo, hi, free):
    """least_squares over the free subset of parameters; returns (v, result)."""
    from scipy.optimize import least_squares
    v = v0.copy()
    idx = np.where(free)[0]

    def resid(z):
        v[idx] = z
        return model(x, v) - y

    z0 = np.clip(v0[idx], lo[idx], hi[idx])
    res = least_squares(resid, z0, bounds=(lo[idx], hi[idx]), method="trf", x_scale="jac",
                        max_nfev=50 + 25 * len(idx))
    v[idx] = res.x
    return v, res


def _covariance(res, m):
    """Covariance of the free parameters, scaled by the reduced chi-square (as lmfit does)."""
    J = res.jac
    n = J.shape[1]
    dof = max(m - n, 1)
    try:
        cov = np.linalg.inv(J.T @ J) * (2.0 * res.cost / dof)
    except np.linalg.LinAlgError:
        return None
    if not np.all(np.isfinite(np.diag(cov))) or np.any(np.diag(cov) < 0):
        return None
    return cov


def fit_group(two_theta, y_sub, group: List[Candidate], wl, doublet=True,
              reach=3.0, min_half_window=0.3, group_id=0, noise=1.0,
              min_r2=C.DEFAULT_MIN_R2, asymmetry=True, min_fwhm=None, log=None) -> List[PeakFit]:
    """
    Fit a group of overlapping peaks. Candidates that turn out to be spurious
    (amplitude not significant, or narrower than ``min_fwhm``) are removed and
    the group is refitted, so that noise cannot distort the real peaks.
    """
    step = float(np.median(np.diff(two_theta)))
    lo_w = min(c.pos - max(reach * c.fwhm_est, min_half_window) for c in group)
    hi_w = max(c.pos + max(reach * c.fwhm_est, min_half_window)
               + (float(ka2_position(c.pos, wl) - c.pos) if doublet else 0.0) for c in group)
    m = (two_theta >= lo_w) & (two_theta <= hi_w)
    x, y = two_theta[m], y_sub[m]
    sw = np.broadcast_to(np.asarray(noise, float), y_sub.shape)[m]
    r = wl["ka2"] / wl["ka1"]
    use_d = bool(doublet and wl.get("ratio", 0) > 0 and r != 1.0)
    ratio = float(wl.get("ratio", 0.0)) if use_d else 0.0

    dropped: List[PeakFit] = []
    group = list(group)
    while True:
        fits, spurious = _fit_once(x, y, sw, group, step, use_d, r, ratio, group_id, min_r2, asymmetry,
                                   min_fwhm)
        if not spurious or len(group) == 1:
            break
        # drop the least significant spurious candidate and refit the rest
        k = min(spurious, key=lambda i: fits[i].snr_fit)
        f = fits[k]
        f.accepted = False
        f.message = f.message if f.message != "ok" else "spurious"
        dropped.append(f)
        if log is not None:
            log.append(f"  {f.center:8.3f}°: removed from its group ({f.message}); group refitted")
        del group[k]
    return fits + dropped


def _fit_once(x, y, sw, group, step, use_d, r, ratio, group_id, min_r2, asymmetry, min_fwhm):
    n = len(group)
    npar = _NP * n + 2
    if len(x) < max(12, npar + 5):
        return [PeakFit(center=c.pos, fwhm=c.fwhm_est, eta=0.5, amplitude=0, height=c.height,
                        accepted=False, flags=["too_few_points"],
                        message=f"too few points ({len(x)}) in window") for c in group], []
    model = _GroupModel(n, use_d, r, ratio)
    v0 = np.zeros(npar)
    lo = np.full(npar, -np.inf)
    hi = np.full(npar, np.inf)
    free = np.ones(npar, bool)
    for k, c in enumerate(group):
        i = 2 + _NP * k
        w = c.fwhm_est
        area0 = max(c.height, 1e-9) * w * (1.2 / (1 + ratio))
        v0[i:i + _NP] = [area0, c.pos, w / 2.0, 0.5, 0.0]
        # narrowest physically possible width (instrument limit); spurious narrow
        # features then finish on this bound and are pruned
        s_min = max(0.5 * float(min_fwhm(c.pos)) if min_fwhm is not None else 0.0, step * 0.25, 1e-4)
        v0[i + 2] = max(v0[i + 2], 1.05 * s_min)
        lo[i:i + _NP] = [0.0, c.pos - (0.6 * w + 2 * step), s_min, 0.0, -0.7]
        hi[i:i + _NP] = [np.inf, c.pos + (0.6 * w + 2 * step), max(6.0 * w, 5 * step, 2.0 * s_min), 1.0, 0.7]
        if not asymmetry:
            lo[i + 4], hi[i + 4], free[i + 4] = -1e-12, 1e-12, False

    try:
        v, res = _solve(model, x, y, v0, lo, hi, free)
    except Exception as e:  # pragma: no cover - defensive
        return [PeakFit(center=c.pos, fwhm=c.fwhm_est, eta=0.5, amplitude=0, height=c.height,
                        accepted=False, flags=["fit_error"], message=str(e)) for c in group], []

    # Parameters that finished on a bound carry no curvature information:
    # fix them and take the covariance of the rest (as the 4.1.1 lmfit retry did).
    at_lo = np.isclose(v, lo, atol=1e-6) & free
    at_hi = np.isclose(v, hi, atol=1e-6) & free
    shape_at_bound = np.zeros(npar, bool)
    for k in range(n):
        i = 2 + _NP * k
        shape_at_bound[i + 3] = at_lo[i + 3] or at_hi[i + 3]
        shape_at_bound[i + 4] = at_lo[i + 4] or at_hi[i + 4]
    free2 = free & ~shape_at_bound
    cov_full = None
    if np.any(shape_at_bound):
        # one refit with the bounded shape parameters fixed gives a clean Jacobian
        try:
            v2, res2 = _solve(model, x, y, v, lo, hi, free2)
            if res2.cost <= res.cost * (1 + 1e-6):
                v, res = v2, res2
            else:
                free2 = free
        except Exception:  # pragma: no cover
            free2 = free
    cov = _covariance(res, len(x))
    names = np.where(free2)[0] if cov is not None and cov.shape[0] == free2.sum() else np.where(free)[0]
    if cov is not None and cov.shape[0] != len(names):
        cov = None
    if cov is not None:
        cov_full = np.full((npar, npar), np.nan)
        cov_full[np.ix_(names, names)] = cov

    resid = model(x, v) - y
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    nfree = int(free2.sum())
    chi2r = float(np.sum((resid / np.maximum(sw, 1e-12)) ** 2) / max(len(x) - nfree, 1))
    xd = np.linspace(x[0], x[-1], 600)
    total = model(xd, v)
    bgd = v[0] + v[1] * xd
    comps = model.components(xd, v)
    success = bool(res.status > 0)

    fits, spurious = [], []
    centers = [v[2 + _NP * k + 1] for k in range(n)]
    for k, c in enumerate(group):
        i = 2 + _NP * k
        amp, cen, sigma, eta, skew = v[i:i + _NP]
        fwhm = 2.0 * sigma
        p1, p2 = comps[k]
        height = amp / pv_integral_breadth(fwhm, eta) if fwhm > 0 else 0.0
        f = PeakFit(center=float(cen), fwhm=float(fwhm), eta=float(eta), amplitude=float(amp),
                    height=float(height), skew=float(skew), r2=r2, chi2_red=chi2r, group_id=group_id, group_size=n,
                    doublet=use_d, data_x=x, data_y=y, fit_x=xd, fit_total=total, fit_peak=p1 + p2,
                    fit_peak_ka1=p1, fit_bg=bgd, residual=resid, snr=c.snr)
        if cov_full is not None and np.isfinite(cov_full[i + 2, i + 2]):
            var_s = max(cov_full[i + 2, i + 2], 0.0)
            f.fwhm_err = float(2 * np.sqrt(var_s))
            if np.isfinite(cov_full[i + 1, i + 1]):
                f.center_err = float(np.sqrt(max(cov_full[i + 1, i + 1], 0.0)))
            if np.isfinite(cov_full[i + 3, i + 3]):
                var_e, cov_se = max(cov_full[i + 3, i + 3], 0.0), cov_full[i + 2, i + 3]
            else:                      # eta fixed at a bound
                var_e, cov_se = 0.05 ** 2, 0.0
            f.eta_err = float(np.sqrt(var_e))
            f.cov_fwhm_eta = np.array([[4 * var_s, 2 * cov_se], [2 * cov_se, var_e]])
            var_a = cov_full[i, i]
            f.snr_fit = float(amp / np.sqrt(var_a)) if np.isfinite(var_a) and var_a > 0 else np.inf
        else:
            f.flags.append("no_uncertainty")
            f.snr_fit = np.nan
        # quality flags
        if eta <= 0.005 or eta >= 0.995:
            f.flags.append("eta_at_bound")
        if sigma >= 0.99 * hi[i + 2] or sigma <= 1.01 * lo[i + 2]:
            f.flags.append("width_at_bound")
        others = [abs(cen - cc) for j, cc in enumerate(centers) if j != k]
        if others and min(others) < 1.5 * fwhm:
            f.flags.append("overlapped")
        if not success:
            f.flags.append("not_converged")
        narrow = None
        if min_fwhm is not None:
            mf = float(min_fwhm(cen))
            if fwhm < mf * 1.01:
                narrow = mf
                f.flags.append("narrower_than_instrument")
        # acceptance
        reasons = []
        # Fit adequacy: R² is low for weak peaks even when the fit is perfect, so a fit
        # whose residuals are consistent with counting noise is also adequate.
        noise_ok = chi2r <= C.MAX_CHI2_RED_WEAK
        if r2 < min_r2 and not noise_ok:
            reasons.append(f"R² = {r2:.3f} < {min_r2} and residuals exceed counting noise (χ²ν = {chi2r:.2f})")
        if "width_at_bound" in f.flags:
            reasons.append("width hit its bound")
        if amp <= 0:
            reasons.append("zero amplitude")
        if not (x[0] < cen < x[-1]):
            reasons.append("centre outside window")
        if "not_converged" in f.flags:
            reasons.append("fit did not converge")
        if np.isfinite(f.snr_fit) and f.snr_fit < C.MIN_AREA_SIGNIFICANCE_ACCEPT:
            reasons.append(f"too weak for a reliable width (area {f.snr_fit:.1f}σ < {C.MIN_AREA_SIGNIFICANCE_ACCEPT:g}σ)")
        if narrow is not None:
            reasons.append(f"narrower than the instrument allows (FWHM {fwhm:.4f}° < {narrow:.4f}°): noise or a spike")
        f.accepted = not reasons
        f.message = "; ".join(reasons) if reasons else "ok"
        # candidates to prune from the group before refitting
        if (narrow is not None or amp <= 0 or "width_at_bound" in f.flags
                or (np.isfinite(f.snr_fit) and f.snr_fit < C.MIN_AMPLITUDE_SIGNIFICANCE)):
            spurious.append(k)
        fits.append(f)
    return fits, spurious


def fit_all(two_theta, y_sub, cands, wl, doublet=True, min_r2=C.DEFAULT_MIN_R2,
            noise=1.0, log=None, asymmetry=True, min_fwhm=None) -> List[PeakFit]:
    groups = group_candidates(cands, wl, doublet=doublet)
    fits = []
    for gid, g in enumerate(groups):
        fits.extend(fit_group(two_theta, y_sub, g, wl, doublet=doublet, group_id=gid,
                              noise=noise, min_r2=min_r2, asymmetry=asymmetry, min_fwhm=min_fwhm,
                              log=log))
    # remove duplicates (two candidates converging on the same peak)
    fits.sort(key=lambda f: f.center)
    out = []
    for f in fits:
        if out and abs(f.center - out[-1].center) < 0.25 * min(f.fwhm, out[-1].fwhm):
            keep = f if f.r2 > out[-1].r2 else out[-1]
            if log is not None:
                log.append(f"  {f.center:8.3f}°: duplicate detection merged")
            out[-1] = keep
        else:
            out.append(f)
    return out
