import numpy as np
from lmfit.models import PseudoVoigtModel
from scipy.special import voigt_profile

from volva import profiles as P


def _num_fwhm(x, y):
    h = y.max() / 2
    above = x[y >= h]
    return above[-1] - above[0]


def test_pseudo_voigt_fwhm_is_two_sigma_for_all_eta():
    """Völva 4.0 used 2.3548σ for the Gaussian part (up to +17.7% error)."""
    m = PseudoVoigtModel()
    x = np.linspace(-3, 3, 600001)
    for eta in (0, 0.3, 0.7, 1):
        y = m.eval(x=x, amplitude=1, center=0, sigma=0.1, fraction=eta)
        assert abs(_num_fwhm(x, y) - P.pv_fwhm(0.1)) < 1e-4


def test_pv_integral_breadth_matches_numeric():
    m = PseudoVoigtModel()
    x = np.linspace(-200, 200, 4000001)
    for eta in (0, 0.5, 1):
        y = m.eval(x=x, amplitude=1, center=0, sigma=0.1, fraction=eta)
        beta = getattr(np, "trapezoid", getattr(np, "trapz", None))(y, x) / y.max()
        assert abs(P.pv_integral_breadth(0.2, eta) - beta) / beta < 2e-3


def test_tch_roundtrip():
    for G, L in [(0.1, 0.0), (0.0, 0.1), (0.05, 0.05), (0.02, 0.3), (0.3, 0.02)]:
        f, eta = P.compose(G, L)
        G2, L2 = P.decompose(f, eta)
        assert abs(G2 - G) < 2e-3 * max(f, 1e-9) + 1e-9
        assert abs(L2 - L) < 2e-3 * max(f, 1e-9) + 1e-9


def test_voigt_integral_breadth_matches_numeric():
    x = np.linspace(-400, 400, 4000001)
    for G, L in [(0.2, 0.0), (0.0, 0.2), (0.1, 0.1), (0.05, 0.3)]:
        sig = G / (2 * np.sqrt(2 * np.log(2)))
        gam = L / 2
        y = voigt_profile(x, max(sig, 1e-12), max(gam, 1e-12))
        beta = 1.0 / y.max()   # unit area
        assert abs(P.voigt_integral_breadth(G, L) - beta) / beta < 2e-3


def test_correction_rules():
    Gs, Ls = P.correct_components(0.5, 0.4, 0.3, 0.1)
    assert abs(Gs - 0.4) < 1e-12 and abs(Ls - 0.3) < 1e-12
    Gs, Ls = P.correct_components(0.1, 0.05, 0.2, 0.1)
    assert Gs == 0 and Ls == 0
