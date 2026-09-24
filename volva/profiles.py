"""
Line-profile mathematics.

Conventions
-----------
* lmfit ``PseudoVoigtModel``: both components share FWHM = 2*sigma; eta is the
  Lorentzian fraction (``fraction``).
* G, L: Gaussian and Lorentzian FWHM of the equivalent Voigt function, found by
  inverting the Thompson-Cox-Hastings (TCH) relations
  (Thompson, Cox & Hastings, J. Appl. Cryst. 20, 79-83, 1987).
* Integral breadth beta = area / height.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq
from scipy.special import erfcx

LN2 = np.log(2.0)
_BG_FACTOR = 0.5 * np.sqrt(np.pi / LN2)   # Gaussian: beta = 1.0645 * FWHM
_BL_FACTOR = 0.5 * np.pi                  # Lorentzian: beta = 1.5708 * FWHM


# ---------------------------------------------------------------------------
# pseudo-Voigt (lmfit parameterisation)
# ---------------------------------------------------------------------------
def pv_fwhm(sigma):
    """FWHM of lmfit's PseudoVoigtModel: 2*sigma for any eta."""
    return 2.0 * np.asarray(sigma)


def pv_integral_breadth(fwhm, eta):
    """Integral breadth of a pseudo-Voigt with given FWHM and Lorentzian fraction."""
    fwhm = np.asarray(fwhm, float)
    eta = np.asarray(eta, float)
    h_per_area = (1 - eta) * 2.0 * np.sqrt(LN2 / np.pi) / fwhm + eta * 2.0 / (np.pi * fwhm)
    return 1.0 / h_per_area


# ---------------------------------------------------------------------------
# TCH relations and their inversion
# ---------------------------------------------------------------------------
_TCH = (1.0, 2.69269, 2.42843, 4.47163, 0.07842, 1.0)


def tch_fwhm(G, L):
    """Total FWHM of a Voigt approximated by TCH, from Gaussian G and Lorentzian L FWHMs."""
    G = np.asarray(G, float)
    L = np.asarray(L, float)
    c = _TCH
    s = (c[0] * G**5 + c[1] * G**4 * L + c[2] * G**3 * L**2 + c[3] * G**2 * L**3
         + c[4] * G * L**4 + c[5] * L**5)
    return np.power(np.maximum(s, 0.0), 0.2)


def tch_eta(q):
    """TCH Lorentzian fraction eta for q = L / FWHM."""
    q = np.asarray(q, float)
    return 1.36603 * q - 0.47719 * q**2 + 0.11116 * q**3


def _build_tables(n=4001):
    q = np.linspace(0.0, 1.0, n)
    eta = tch_eta(q)                      # monotonic on [0, 1], eta(1) = 1
    r = np.empty_like(q)                  # r = G / FWHM for each q
    for i, qi in enumerate(q):
        if qi >= 1.0:
            r[i] = 0.0
        else:
            r[i] = brentq(lambda g: tch_fwhm(g, qi) - 1.0, 0.0, 1.0 + 1e-12)
    return q, eta, r


_Q, _ETA, _R = _build_tables()


def decompose(fwhm, eta):
    """
    pseudo-Voigt (FWHM, eta) -> (G, L) Voigt component FWHMs. Vectorised.

    eta is clipped to [0, 1]. Returns arrays with the broadcast shape.
    """
    fwhm = np.asarray(fwhm, float)
    eta = np.clip(np.asarray(eta, float), 0.0, 1.0)
    q = np.interp(eta, _ETA, _Q)
    r = np.interp(q, _Q, _R)
    return r * fwhm, q * fwhm


def compose(G, L):
    """(G, L) -> (FWHM, eta) using TCH."""
    f = tch_fwhm(G, L)
    with np.errstate(invalid="ignore", divide="ignore"):
        q = np.where(f > 0, np.asarray(L, float) / f, 0.0)
    return f, tch_eta(q)


# ---------------------------------------------------------------------------
# Voigt integral breadth
# ---------------------------------------------------------------------------
def voigt_integral_breadth(G, L):
    """
    Integral breadth of a Voigt from its Gaussian and Lorentzian FWHMs.

    beta = beta_G / erfcx(k), k = beta_L / (sqrt(pi) * beta_G)
    (e.g. Langford, J. Appl. Cryst. 11, 10-14, 1978).
    """
    G = np.asarray(G, float)
    L = np.asarray(L, float)
    bG = _BG_FACTOR * G
    bL = _BL_FACTOR * L
    with np.errstate(divide="ignore", invalid="ignore"):
        k = np.where(bG > 0, bL / (np.sqrt(np.pi) * np.where(bG > 0, bG, 1.0)), np.inf)
        beta = np.where(bG > 0, bG / erfcx(np.where(np.isfinite(k), k, 0.0)), bL)
    return beta


def gaussian_breadth(G):
    return _BG_FACTOR * np.asarray(G, float)


def lorentzian_breadth(L):
    return _BL_FACTOR * np.asarray(L, float)


def voigt_fwhm(G, L):
    """Olivero & Longbothum (1977) Voigt FWHM, used only for display."""
    G = np.asarray(G, float)
    L = np.asarray(L, float)
    return 0.5346 * L + np.sqrt(0.2166 * L**2 + G**2)


# ---------------------------------------------------------------------------
# Instrumental correction of component widths
# ---------------------------------------------------------------------------
def correct_components(G_obs, L_obs, G_inst, L_inst):
    """
    Remove the instrument profile from an observed Voigt.

    Voigt convolution: Gaussian FWHMs add in quadrature, Lorentzian FWHMs add
    linearly. Negative results are clipped to zero (not resolvable).
    """
    G_s = np.sqrt(np.clip(np.asarray(G_obs) ** 2 - np.asarray(G_inst) ** 2, 0.0, None))
    L_s = np.clip(np.asarray(L_obs) - np.asarray(L_inst), 0.0, None)
    return G_s, L_s
