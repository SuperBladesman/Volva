"""Prototype: Voigt-aware instrumental correction (TCH decomposition), for audit comparison only."""
import numpy as np
from scipy.optimize import brentq
def tch_gamma(G, L):
    return (G**5 + 2.69269*G**4*L + 2.42843*G**3*L**2 + 4.47163*G**2*L**3 + 0.07842*G*L**4 + L**5)**0.2
def decompose(fwhm, eta):
    """pseudo-Voigt (FWHM, eta) -> (Gaussian FWHM, Lorentzian FWHM) via TCH (Thompson, Cox & Hastings 1987)."""
    eta = min(max(eta, 0.0), 1.0)
    f = lambda q: 1.36603*q - 0.47719*q**2 + 0.11116*q**3 - eta
    q = 0.0 if eta <= 0 else (1.0 if eta >= 1 else brentq(f, 0, 1))
    L = q*fwhm
    G = brentq(lambda g: tch_gamma(g, L) - fwhm, 0, fwhm*1.0001) if q < 1 else 0.0
    return G, L
