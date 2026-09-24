"""Physical constants, reference data and default settings."""

from __future__ import annotations

# ---------------------------------------------------------------------------
# X-ray wavelengths (nm). International Tables for Crystallography Vol. C,
# Table 4.2.2.1. Ratio = I(Ka2)/I(Ka1).
# Wavelengths read from an XRDML file override these values.
# ---------------------------------------------------------------------------
WAVELENGTHS = {
    "Cu Ka": {"ka1": 0.1540562, "ka2": 0.1544390, "ratio": 0.5},
    "Co Ka": {"ka1": 0.1788965, "ka2": 0.1792850, "ratio": 0.5},
    "Fe Ka": {"ka1": 0.1936042, "ka2": 0.1939980, "ratio": 0.5},
    "Cr Ka": {"ka1": 0.2289700, "ka2": 0.2293606, "ratio": 0.5},
    "Mo Ka": {"ka1": 0.0709300, "ka2": 0.0713590, "ratio": 0.5},
    "Ag Ka": {"ka1": 0.0559407, "ka2": 0.0563798, "ratio": 0.5},
}

# ---------------------------------------------------------------------------
# Size conventions. Völva measures INTEGRAL BREADTH (area / height) of the
# sample profile. With K = 1 the Scherrer equation then gives the
# volume-weighted mean column length L_vol (Stokes & Wilson, 1942), which is
# the quantity certified for NIST SRM 1979. For spheres of diameter D,
# L_vol = 3D/4, so D = (4/3) L_vol.
# ---------------------------------------------------------------------------
SIZE_CONVENTIONS = {
    "L_vol (volume-weighted column length, K = 1)": 1.0,
    "Sphere diameter (K = 4/3)": 4.0 / 3.0,
}
DEFAULT_SIZE_CONVENTION = "L_vol (volume-weighted column length, K = 1)"

# ---------------------------------------------------------------------------
# Line-profile standards (cubic). Used only for a position sanity check
# (wavelength / zero-shift), never for the width calibration itself.
# ---------------------------------------------------------------------------
STANDARDS = {
    "LaB6 (NIST SRM 660c)": {"a_nm": 0.4156826, "lattice": "cubic_P"},
    "Si (NIST SRM 640 series)": {"a_nm": 0.5431, "lattice": "cubic_Fd3m"},
    "Other / skip position check": None,
}
DEFAULT_STANDARD = "LaB6 (NIST SRM 660c)"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
DEFAULT_MIN_R2 = 0.90            # minimum R^2 for an accepted sample peak fit
DEFAULT_STD_MIN_R2 = 0.95        # minimum R^2 for an accepted standard peak
DEFAULT_PROMINENCE_SIGMA = 5.0   # detection: net peak area above a local baseline, in counting-noise σ
MIN_AMPLITUDE_SIGNIFICANCE = 3.0  # below this a group member is treated as spurious and removed
MIN_AREA_SIGNIFICANCE_ACCEPT = 5.0  # an accepted peak needs its fitted area >= this many standard errors
MIN_FWHM_FRACTION_OF_INSTRUMENT = 0.7  # sample peaks narrower than this x instrument FWHM are spurious
MAX_CHI2_RED_WEAK = 1.5         # a fit with residuals within counting noise is adequate even if R² is low
DEFAULT_MIN_REL_HEIGHT = 0.01    # detection floor, fraction of the strongest peak
DEFAULT_RESOLUTION_RATIO = 1.2   # beta_obs / beta_inst below this -> resolution-limited
DEFAULT_MC_SAMPLES = 2000        # Monte Carlo draws for uncertainty propagation
DEFAULT_SEED = 20260924          # fixed seed: identical input -> identical output
MIN_PEAKS_LINEAR = 3             # minimum peaks for W-H, H-W, Voigt W-H
MIN_PEAKS_MONSHI = 3

# Background (SNIP) half-width limits, degrees 2theta
BG_MIN_HALF_WIDTH = 0.8
BG_MAX_HALF_WIDTH = 6.0
BG_WIDTH_FACTOR = 2.5            # half-width = factor x widest typical FWHM

# Default instrument profile uncertainty floor (relative) used if the
# calibration residuals are smaller than this.
INST_REL_UNC_FLOOR = 0.03
ETA_BOUND_UNC = 0.15        # uncertainty assigned to a standard peak whose eta hit 0 or 1 (censored)
