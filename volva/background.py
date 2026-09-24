"""
Background estimation: SNIP (Ryan et al., Nucl. Instrum. Methods B34, 396, 1988)
with the clipping window defined in DEGREES 2θ, so results do not depend on
the scan step size.
"""

from __future__ import annotations

import numpy as np
from scipy.signal import savgol_filter

from . import constants as C


def _lls(y):
    return np.log(np.log(np.sqrt(np.maximum(y, 0) + 1.0) + 1.0) + 1.0)


def _lls_inv(z):
    return (np.exp(np.exp(z) - 1.0) - 1.0) ** 2 - 1.0


def snip(two_theta: np.ndarray, intensity: np.ndarray, half_width_deg: float) -> np.ndarray:
    """
    SNIP background with a clipping half-width given in degrees.

    Works in LLS space (log-log-sqrt), which follows Poisson statistics and
    avoids clipping into the flanks of intense peaks. Vectorised.
    """
    step = float(np.median(np.diff(two_theta)))
    m = max(1, int(round(half_width_deg / step)))
    n = len(intensity)
    m = min(m, (n - 1) // 2)
    z = _lls(intensity)
    for p in range(1, m + 1):
        mid = 0.5 * (z[:-2 * p] + z[2 * p:])
        z[p:n - p] = np.minimum(z[p:n - p], mid)
    bg = _lls_inv(z)
    return np.minimum(bg, intensity)


def noise_sigma(intensity: np.ndarray) -> float:
    """Robust noise estimate from first differences (MAD), in intensity units."""
    d = np.diff(intensity)
    mad = np.median(np.abs(d - np.median(d)))
    s = 1.4826 * mad / np.sqrt(2.0)
    return float(s) if s > 0 else float(np.std(d) / np.sqrt(2.0) or 1.0)


def local_noise(intensity: np.ndarray, background: np.ndarray) -> np.ndarray:
    """
    Point-by-point noise standard deviation, sqrt(gain * expected counts).

    Counting noise grows with intensity, so a single global sigma under-states
    the noise where the background is high (e.g. at low angle) and over-states
    it where it is low. The gain is estimated from the data (1 for raw counts),
    so scaled or cps data are handled too.
    """
    I = np.asarray(intensity, float)
    z = np.diff(I) / np.sqrt(np.maximum(np.abs(I[:-1]) + np.abs(I[1:]), 1e-12))
    gain = (1.4826 * np.median(np.abs(z - np.median(z)))) ** 2
    if not np.isfinite(gain) or gain <= 0:
        return np.full_like(I, noise_sigma(I))
    floor = max(float(np.median(np.abs(I))) * 1e-3, 1e-12)
    # expected counts: smoothed total intensity (peaks included), never below the background
    n = len(I)
    w = min(7, (n // 2) * 2 - 1)
    level = savgol_filter(I, w, 2) if w >= 5 else I
    level = np.maximum(level, np.asarray(background, float))
    return np.sqrt(gain * np.maximum(level, floor))


def estimate_fwhms(two_theta, y_sub, n_max=12):
    """Rough FWHM estimates (degrees) of the strongest features in background-subtracted data."""
    from scipy.signal import find_peaks
    step = float(np.median(np.diff(two_theta)))
    w = max(5, int(0.1 / step) | 1)
    w = min(w, (len(y_sub) // 2) * 2 - 1)
    ys = savgol_filter(y_sub, w, 2) if w >= 5 else y_sub
    sig = noise_sigma(y_sub)
    idx, props = find_peaks(ys, prominence=8 * sig, width=2)
    if len(idx) == 0:
        return np.array([])
    order = np.argsort(props["prominences"])[::-1][:n_max]
    return props["widths"][order] * step


def auto_half_width(two_theta, intensity) -> tuple[float, str]:
    """
    Choose the SNIP half-width from the peak widths in the pattern.

    Pass 1: a wide (4°) SNIP to expose peaks. Estimate FWHMs of the strongest
    features. Pass 2 half-width = BG_WIDTH_FACTOR x 75th-percentile FWHM,
    clipped to [BG_MIN_HALF_WIDTH, BG_MAX_HALF_WIDTH].
    """
    bg0 = snip(two_theta, intensity, 4.0)
    fw = estimate_fwhms(two_theta, intensity - bg0)
    if len(fw) == 0:
        return 2.0, "no clear peaks found; using 2.0°"
    typical = float(np.percentile(fw, 75))
    hw = float(np.clip(C.BG_WIDTH_FACTOR * typical, C.BG_MIN_HALF_WIDTH, C.BG_MAX_HALF_WIDTH))
    note = f"auto: 75th-percentile FWHM ≈ {typical:.3f}° → half-width {hw:.2f}°"
    if C.BG_WIDTH_FACTOR * typical > C.BG_MAX_HALF_WIDTH:
        note += " (capped; very broad peaks — check the background plot)"
    return hw, note
