"""Method algebra on exact, noise-free corrected peaks."""
import numpy as np

from volva.correction import CorrectedPeak
from volva.methods import halder_wagner, ols, run_methods, voigt_wh, williamson_hall
from volva.peaks import PeakFit

LAM = 0.1540562


def make_peaks(D, e, tts=(28, 33, 47, 56, 69, 76, 88, 95, 107, 115), lorentz_size=True):
    peaks = []
    for tt in tts:
        th = np.deg2rad(tt / 2)
        bL = LAM / (D * np.cos(th))                 # K = 1, Lorentzian size breadth
        bG = 4 * e * np.tan(th)                     # Gaussian strain breadth
        from volva.profiles import voigt_integral_breadth
        # integral breadth of the Voigt with these component breadths
        G = bG / (0.5 * np.sqrt(np.pi / np.log(2)))
        L = bL / (np.pi / 2)
        b = float(voigt_integral_breadth(G, L))
        f = PeakFit(center=tt, fwhm=0.1, eta=0.5, amplitude=1, height=1, accepted=True)
        p = CorrectedPeak(fit=f, two_theta=tt, theta=th, d_nm=LAM / (2 * np.sin(th)))
        p.beta_s, p.beta_Ls, p.beta_Gs = b, bL, bG
        p.beta_s_err = p.beta_Ls_err = p.beta_Gs_err = 1e-9
        n = 50
        p.mc = {"beta_s": np.full(n, b), "beta_Ls": np.full(n, bL), "beta_Gs": np.full(n, bG)}
        peaks.append(p)
    return peaks


def test_voigt_wh_recovers_exact_inputs():
    r = voigt_wh(make_peaks(40.0, 1e-3), LAM, 1.0)
    assert abs(r.size_nm - 40.0) / 40 < 1e-6
    assert abs(r.strain - 1e-3) / 1e-3 < 1e-6


def test_scherrer_exact_when_strain_free():
    r = run_methods(make_peaks(40.0, 0.0), LAM, 1.0)
    assert abs(r["scherrer"].size_nm - 40.0) < 1e-6
    assert abs(r["single_line"].size_nm - 40.0) < 1e-6
    assert abs(r["monshi"].slope - 1.0) < 1e-6


def test_hw_recovers_voigt_inputs_approximately():
    """H-W uses the parabolic Voigt approximation; expect a few % agreement."""
    r = halder_wagner(make_peaks(40.0, 1e-3), LAM, 1.0)
    assert abs(r.size_nm - 40.0) / 40 < 0.06
    assert abs(r.strain - 1e-3) / 1e-3 < 0.15


def test_ssp_dimensionally_consistent_form_equals_hw():
    """(d β cosθ / λ)² = (K/D)(d² β cosθ / λ) + (2e)²  is the H-W regression rescaled (Pal 2025)."""
    peaks = make_peaks(40.0, 1e-3)
    hw = halder_wagner(peaks, LAM, 1.0)
    x = np.array([p.d_nm**2 * p.beta_s * np.cos(p.theta) / LAM for p in peaks])
    y = np.array([(p.d_nm * p.beta_s * np.cos(p.theta) / LAM) ** 2 for p in peaks])
    m, c, *_ = ols(x, y)
    assert abs(1.0 / m - hw.size_nm) / hw.size_nm < 1e-9
    assert abs(np.sqrt(c) / 2 - hw.strain) / hw.strain < 1e-9


def test_wh_bias_direction_with_gaussian_strain():
    """Voigt breadths are smaller than the linear sum, so UDM W-H gives size too large, strain too small."""
    r = williamson_hall(make_peaks(40.0, 1e-3), LAM, 1.0)
    assert r.size_nm > 40.0      # linear-sum W-H on Voigt breadths: intercept too small -> size too big
    assert r.strain < 1e-3


def test_negative_slope_not_reported_as_strain():
    peaks = make_peaks(40.0, 0.0)
    for i, p in enumerate(peaks):         # breadth falling with angle
        p.beta_s *= 1.3 - 0.05 * i
    r = williamson_hall(peaks, LAM, 1.0)
    assert r.strain_significant is False
    assert any("Negative slope" in w for w in r.warnings)
