"""
Instrumental correction of fitted sample peaks, with Monte Carlo propagation
of fit and calibration uncertainties.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
from scipy.optimize import brentq

from .instrument import InstrumentProfile
from .peaks import PeakFit
from .profiles import (correct_components, decompose, gaussian_breadth,
                       lorentzian_breadth, voigt_integral_breadth)

DEG = np.pi / 180.0


@dataclass
class CorrectedPeak:
    fit: PeakFit
    two_theta: float
    theta: float                   # rad
    d_nm: float
    hkl: str = ""
    use: bool = True
    # component FWHMs (degrees 2theta)
    G_obs: float = 0.0
    L_obs: float = 0.0
    G_inst: float = 0.0
    L_inst: float = 0.0
    G_s: float = 0.0
    L_s: float = 0.0
    # integral breadths
    beta_obs_deg: float = 0.0
    beta_inst_deg: float = 0.0
    beta_s: float = 0.0            # rad, total sample (Voigt)
    beta_Ls: float = 0.0           # rad, Lorentzian (size-like) part
    beta_Gs: float = 0.0           # rad, Gaussian (strain-like) part
    beta_s_err: float = 0.0
    beta_Ls_err: float = 0.0
    beta_Gs_err: float = 0.0
    ratio: float = np.inf          # beta_obs / beta_inst
    resolution_limited: bool = False
    beta_min: float = np.nan       # rad: sample breadth at the resolution threshold (instrument property)
    beta_upper: float = np.nan     # rad: 2-sigma upper limit on the sample breadth (data)
    component_mode: str = "both"
    flags: List[str] = field(default_factory=list)
    mc: dict = field(default_factory=dict, repr=False)

    @property
    def lorentz_fraction(self) -> float:
        s = self.beta_Ls + self.beta_Gs
        return self.beta_Ls / s if s > 0 else np.nan


def _min_resolvable_breadth(G_i, L_i, ratio):
    """Lorentzian sample breadth (deg) that makes beta_obs = ratio * beta_inst."""
    ratio = max(float(ratio), 1.01)
    b_inst = float(voigt_integral_breadth(G_i, L_i))
    target = ratio * b_inst
    f = lambda Ladd: float(voigt_integral_breadth(G_i, L_i + Ladd)) - target  # noqa: E731
    hi = max(b_inst, 1e-4)
    while f(hi) < 0:
        hi *= 2
    L_add = brentq(f, 0.0, hi)
    return float(lorentzian_breadth(L_add))


def _bisect(fun, target, hi, n_iter=48):
    """Vectorised bisection for monotonic increasing fun(x) = target on [0, hi]."""
    lo = np.zeros_like(target)
    hi = np.array(hi, float) + np.zeros_like(target)
    for _ in range(8):                      # make sure the bracket is wide enough
        grow = fun(hi) < target
        if not np.any(grow):
            break
        hi = np.where(grow, hi * 4, hi)
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        below = fun(mid) < target
        lo = np.where(below, mid, lo)
        hi = np.where(below, hi, mid)
    return 0.5 * (lo + hi)


def constrained_correction(G_o, L_o, G_i, L_i, only=None):
    """
    Sample (G_s, L_s) from observed and instrument Voigt components, preserving
    the observed TOTAL integral breadth when one component is smaller than the
    instrument's.

    Noise in the TCH split of an observed profile is anti-correlated between G
    and L. Clipping each difference at zero separately (Völva 4.1 pre-release)
    kept the excess in one component while discarding the deficit in the other,
    inflating the sample breadth for large crystallites. Here, if G_o < G_i the
    sample is taken as purely Lorentzian with the width that reproduces the
    observed breadth, and vice versa.
    """
    scalar = all(np.ndim(v) == 0 for v in (G_o, L_o, G_i, L_i))
    G_o, L_o, G_i, L_i = (np.array(v, float) for v in np.broadcast_arrays(
        *(np.atleast_1d(np.asarray(v, float)) for v in (G_o, L_o, G_i, L_i))))
    G_s = np.atleast_1d(np.sqrt(np.clip(G_o**2 - G_i**2, 0.0, None)))
    L_s = np.atleast_1d(np.clip(L_o - L_i, 0.0, None))
    b_obs = voigt_integral_breadth(G_o, L_o)
    b_inst = voigt_integral_breadth(G_i, L_i)
    excess = b_obs > b_inst
    if only == "L":          # Gaussian excess not significant: whole excess breadth is Lorentzian
        gdef = np.ones_like(G_o, bool)
        ldef = np.zeros_like(G_o, bool)
        both = np.zeros_like(G_o, bool)
    elif only == "G":        # Lorentzian excess not significant: whole excess breadth is Gaussian
        ldef = np.ones_like(G_o, bool)
        gdef = np.zeros_like(G_o, bool)
        both = np.zeros_like(G_o, bool)
    else:
        gdef = (G_o < G_i) & (L_o >= L_i)
        ldef = (L_o < L_i) & (G_o >= G_i)
        both = (G_o < G_i) & (L_o < L_i)
    if np.any(gdef):
        m = gdef & excess
        if np.any(m):
            Gi_m, Li_m, tgt = G_i[m], L_i[m], b_obs[m]
            L_s[m] = _bisect(lambda x: voigt_integral_breadth(Gi_m, Li_m + x), tgt, np.maximum(tgt, 1e-6))
        G_s[gdef] = 0.0
        L_s[gdef & ~excess] = 0.0
    if np.any(ldef):
        m = ldef & excess
        if np.any(m):
            Gi_m, Li_m, tgt = G_i[m], L_i[m], b_obs[m]
            G_s[m] = _bisect(lambda x: voigt_integral_breadth(np.sqrt(Gi_m**2 + x**2), Li_m), tgt,
                             np.maximum(tgt, 1e-6))
        L_s[ldef] = 0.0
        G_s[ldef & ~excess] = 0.0
    G_s[both] = 0.0
    L_s[both] = 0.0
    # Both components above the instrument: rescale the (G_s, L_s) pair so the
    # reconstructed profile reproduces the observed total breadth. In the
    # noise-free case the scale is exactly 1; with noise it stops the quadrature
    # term sqrt(G_o^2 - G_i^2) from amplifying a tiny excess into a large breadth.
    reg = ~(gdef | ldef | both)
    if np.any(reg):
        Gs0, Ls0, Gi_r, Li_r, tgt = G_s[reg], L_s[reg], G_i[reg], L_i[reg], b_obs[reg]
        ok = tgt > voigt_integral_breadth(Gi_r, Li_r)
        scale = np.zeros_like(tgt)
        if np.any(ok):
            fun = lambda s: voigt_integral_breadth(np.sqrt(Gi_r[ok]**2 + (s * Gs0[ok])**2),  # noqa: E731
                                                   Li_r[ok] + s * Ls0[ok])
            scale[ok] = _bisect(fun, tgt[ok], np.full(ok.sum(), 2.0))
        G_s[reg] = scale * Gs0
        L_s[reg] = scale * Ls0
    if scalar:
        return float(G_s[0]), float(L_s[0])
    return G_s, L_s


class _ShiftedProfile:
    """
    Instrument profile shifted by k standard deviations (systematic sensitivity).

    mode "offset": the same shift at every angle.
    mode "tilt":   shift varying linearly with tan(theta) from -k to +k across the
                   calibrated range, which mimics (or hides) strain.
    """

    def __init__(self, prof, k, mode="offset"):
        self.p, self.k, self.mode = prof, k, mode
        lo, hi = prof.two_theta_range
        self.t0 = np.tan(np.deg2rad(lo / 2))
        self.t1 = np.tan(np.deg2rad(hi / 2))

    def _k(self, tt):
        if self.mode == "offset" or self.t1 <= self.t0:
            return self.k
        t = np.tan(np.deg2rad(np.asarray(tt, float) / 2))
        return self.k * np.clip(2 * (t - self.t0) / (self.t1 - self.t0) - 1, -1.5, 1.5)

    def G_at(self, tt):
        return np.clip(self.p.G_at(tt) + self._k(tt) * self.p.sigma_G(tt), 0, None)

    def L_at(self, tt):
        return np.clip(self.p.L_at(tt) + self._k(tt) * self.p.sigma_L(tt), 0, None)

    def sigma_G(self, tt):
        return 0.0 * np.asarray(tt, float)

    def sigma_L(self, tt):
        return 0.0 * np.asarray(tt, float)

    def covers(self, tt, margin=2.0):
        return self.p.covers(tt, margin)


def correct_peaks(fits: List[PeakFit], wl: dict, profile: Optional[InstrumentProfile],
                  resolution_ratio: float = 1.2, n_mc: int = 2000, seed: int = 0) -> List[CorrectedPeak]:
    rng = np.random.default_rng(seed)
    resolution_ratio = max(float(resolution_ratio), 1.01)
    out = []
    for f in fits:
        if not f.accepted:
            continue
        th = np.deg2rad(f.center / 2)
        cp = CorrectedPeak(fit=f, two_theta=f.center, theta=th,
                           d_nm=wl["ka1"] / (2 * np.sin(th)))
        G_o, L_o = decompose(f.fwhm, f.eta)
        cp.G_obs, cp.L_obs = float(G_o), float(L_o)
        cp.beta_obs_deg = float(voigt_integral_breadth(G_o, L_o))

        if profile is not None:
            G_i = float(profile.G_at(f.center))
            L_i = float(profile.L_at(f.center))
            sG = float(profile.sigma_G(f.center))
            sL = float(profile.sigma_L(f.center))
            if not profile.covers(f.center):
                cp.flags.append("outside_calibrated_range")
        else:
            G_i = L_i = sG = sL = 0.0
            cp.flags.append("no_instrument_correction")
        cp.G_inst, cp.L_inst = G_i, L_i
        cp.beta_inst_deg = float(voigt_integral_breadth(G_i, L_i)) if (G_i + L_i) > 0 else 0.0

        # ---- Monte Carlo draws of the observed components -----------------
        mean = np.array([f.fwhm, f.eta])
        if f.cov_fwhm_eta is not None and np.all(np.isfinite(f.cov_fwhm_eta)):
            cov = f.cov_fwhm_eta.copy() + np.diag([1e-12, 1e-12])
        else:
            cov = np.diag([(0.05 * f.fwhm) ** 2, 0.08 ** 2])
            cp.flags.append("assumed_fit_uncertainty")
        if "eta_at_bound" in f.flags:
            # lmfit's error on a parameter pinned at its bound is meaningless (~1e-5)
            cov[1, 1] = max(cov[1, 1], 0.05 ** 2)
            cov[0, 1] = cov[1, 0] = 0.0
        try:
            s = rng.multivariate_normal(mean, cov, size=n_mc, method="eigh")
        except Exception:
            s = np.column_stack([rng.normal(mean[0], np.sqrt(abs(cov[0, 0])), n_mc),
                                 rng.normal(mean[1], np.sqrt(abs(cov[1, 1])), n_mc)])
        fw = np.clip(s[:, 0], 1e-6, None)
        et = np.clip(s[:, 1], 0.0, 1.0)
        Go, Lo = decompose(fw, et)
        Gi = np.clip(rng.normal(G_i, sG, n_mc), 0, None) if sG > 0 else np.full(n_mc, G_i)
        Li = np.clip(rng.normal(L_i, sL, n_mc), 0, None) if sL > 0 else np.full(n_mc, L_i)

        # A sample component is assigned only if its excess over the instrument is
        # significant (> 2 sigma, fit + calibration). Otherwise the whole excess
        # breadth goes to the other component. This stops the quadrature term
        # sqrt(G_o^2 - G_i^2) turning noise into a large Gaussian breadth.
        sGo, sLo = float(np.std(Go)), float(np.std(Lo))
        g_sig = (G_o - G_i) > 2 * np.hypot(sGo, sG)
        l_sig = (L_o - L_i) > 2 * np.hypot(sLo, sL)
        only = None if (g_sig and l_sig) else ("G" if g_sig else "L")
        if only == "L" and not l_sig:
            only = "L"      # neither significant: treat any excess as size-like (Lorentzian)
        cp.component_mode = {None: "both", "L": "Lorentzian only", "G": "Gaussian only"}[only]
        G_s, L_s = constrained_correction(G_o, L_o, G_i, L_i, only=only)
        cp.G_s, cp.L_s = float(G_s), float(L_s)
        cp.beta_s = float(voigt_integral_breadth(G_s, L_s)) * DEG if (G_s + L_s) > 0 else 0.0
        cp.beta_Ls = float(lorentzian_breadth(L_s)) * DEG
        cp.beta_Gs = float(gaussian_breadth(G_s)) * DEG
        Gs, Ls = constrained_correction(Go, Lo, Gi, Li, only=only)
        b_s = voigt_integral_breadth(Gs, Ls) * DEG
        b_L = lorentzian_breadth(Ls) * DEG
        b_G = gaussian_breadth(Gs) * DEG
        cp.mc = {"beta_s": b_s, "beta_Ls": b_L, "beta_Gs": b_G}
        cp.beta_s_err = float(np.std(b_s))
        cp.beta_Ls_err = float(np.std(b_L))
        cp.beta_Gs_err = float(np.std(b_G))

        # ---- resolution --------------------------------------------------
        # "resolution_limited": sample breadth not statistically resolved
        # (< 2 sigma). Such reflections give only a lower bound on size, taken
        # from the data (beta_s + 2 sigma), and are excluded from means.
        # "near_resolution_limit": beta_obs/beta_inst below the user threshold;
        # kept (excluding them would bias the mean low) but flagged.
        if cp.beta_inst_deg > 0:
            cp.ratio = cp.beta_obs_deg / cp.beta_inst_deg
            cp.beta_min = _min_resolvable_breadth(G_i, L_i, resolution_ratio) * DEG
            if cp.ratio < resolution_ratio:
                cp.flags.append("near_resolution_limit")
        if cp.beta_s < 2 * cp.beta_s_err or cp.beta_s <= 0:
            cp.resolution_limited = True
            cp.flags.append("resolution_limited")
        cp.beta_upper = cp.beta_s + 2 * cp.beta_s_err
        out.append(cp)
    return out
