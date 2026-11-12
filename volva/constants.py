"""Physical constants, reference data and default settings."""

from __future__ import annotations

# ---------------------------------------------------------------------------
# X-ray wavelengths (nm). SOURCED: International Tables for Crystallography
# Vol. C, Table 4.2.2.1. Ratio = I(Ka2)/I(Ka1); 0.5 is the usual approximation,
# not a value from that table.
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
# L_vol = 3D/4, so D = (4/3) L_vol (Langford & Wilson, J. Appl. Cryst. 11,
# 102, 1978). SOURCED.
# ---------------------------------------------------------------------------
SIZE_CONVENTIONS = {
    "L_vol (volume-weighted column length, K = 1)": 1.0,
    "Sphere diameter (K = 4/3)": 4.0 / 3.0,
}
DEFAULT_SIZE_CONVENTION = "L_vol (volume-weighted column length, K = 1)"

# ---------------------------------------------------------------------------
# Line-profile standards (cubic). Used only for a position sanity check
# (wavelength / zero-shift), never for the width calibration itself.
# LaB6: SOURCED, the NIST SRM 660c certificate. Si: a nominal silicon lattice
# parameter, not taken from a particular certificate; adequate for the check,
# which matches peaks within 0.3°.
# ---------------------------------------------------------------------------
STANDARDS = {
    "LaB6 (NIST SRM 660c)": {"a_nm": 0.4156826, "lattice": "cubic_P"},
    "Si (NIST SRM 640 series)": {"a_nm": 0.5431, "lattice": "cubic_Fd3m"},
    "Other / skip position check": None,
}
DEFAULT_STANDARD = "LaB6 (NIST SRM 660c)"

# ---------------------------------------------------------------------------
# Defaults
#
# Each acceptance default is marked SOURCED, with its reference, or UNSOURCED,
# with what it was set on. The README table "Thresholds and where they come
# from" matches this list. Thresholds still written into the code, rather than
# here, are listed in docs/THRESHOLDS.md.
# ---------------------------------------------------------------------------

# UNSOURCED. 0.8 in volva-xrd (3.x defined 0.8 but did not apply it), 0.65 in
# 4.0, 0.90 from 4.1; no measurement or reference is recorded for the value.
# Weak peaks with a lower R² are still accepted if their residuals are within
# counting noise (MAX_CHI2_RED_WEAK).
DEFAULT_MIN_R2 = 0.90            # minimum R^2 for an accepted sample peak fit

# UNSOURCED. Stricter than for samples because standard peaks are strong
# (4.0 used 0.7); no measurement or reference is recorded for the value.
DEFAULT_STD_MIN_R2 = 0.95        # minimum R^2 for an accepted standard peak

# UNSOURCED. 4.0 required a 5σ height-to-noise ratio; 4.1.2 made the test the
# net peak area over the local counting noise and kept 5σ. On four real
# patterns, 3σ and 5σ selected the same peaks (CHANGELOG, 4.1.2).
DEFAULT_PROMINENCE_SIGMA = 5.0   # detection: net peak area above a local baseline, in counting-noise σ

# UNSOURCED (4.1.2). Despite the name, the test is on the fitted area (the
# amplitude of the pseudo-Voigt is its area), in standard errors: a group
# member below it is treated as spurious, removed, and the group refitted.
MIN_AMPLITUDE_SIGNIFICANCE = 3.0  # below this a group member is treated as spurious and removed

# UNSOURCED (4.1.2). Equal to the detection default above; no other reason is
# recorded.
MIN_AREA_SIGNIFICANCE_ACCEPT = 5.0  # an accepted peak needs its fitted area >= this many standard errors

# UNSOURCED (4.1.2). The physical limit is 1: a sample peak cannot be narrower
# than the instrument profile. The margin below 1 allows for the error in a
# fitted width, but no measurement is recorded for 0.7. The limit is never
# less than 1.5 scan steps (pipeline.py).
MIN_FWHM_FRACTION_OF_INSTRUMENT = 0.7  # sample peaks narrower than this x instrument FWHM are spurious

# UNSOURCED (4.1.2). Residuals consistent with counting noise give a reduced
# chi-square near 1; 1.5 is the margin allowed. No measurement is recorded.
MAX_CHI2_RED_WEAK = 1.5         # a fit with residuals within counting noise is adequate even if R² is low

# UNSOURCED. A floor for detection, as a fraction of the strongest peak; the
# 5σ significance test above does the real work. 4.0 detected at 0.05 of the
# maximum (its introduction said 0.01). No measurement is recorded for 0.01.
DEFAULT_MIN_REL_HEIGHT = 0.01    # detection floor, fraction of the strongest peak

# UNSOURCED. Introduced in 3.1 (November 2025) as "sample width < 20 % of
# instrumental" and kept since. A reflection below it is flagged as near the
# resolution limit, and the ratio sets the largest measurable size. If more
# than half the reflections are near the limit or unresolved, the headline
# becomes a 2σ lower bound and no mean is reported (diagnostics.py). Whether a
# single reflection gives a size or only a bound is decided separately, by a
# 2σ test on its sample breadth (correction.py).
DEFAULT_RESOLUTION_RATIO = 1.2   # beta_obs / beta_inst below this -> near the resolution limit

# UNSOURCED. No convergence study is recorded for the number of draws.
DEFAULT_MC_SAMPLES = 2000        # Monte Carlo draws for uncertainty propagation

# Not a threshold: any fixed seed makes the output reproducible.
DEFAULT_SEED = 20260924          # fixed seed: identical input -> identical output

# UNSOURCED. Three points are the fewest that leave one residual degree of
# freedom for a two-parameter straight line.
MIN_PEAKS_LINEAR = 3             # minimum peaks for W-H, H-W, Voigt W-H
MIN_PEAKS_MONSHI = 3

# Background (SNIP) half-width limits, degrees 2theta.
# The SNIP algorithm is SOURCED (Ryan et al., Nucl. Instrum. Methods B 34, 396,
# 1988). These limits and the factor are UNSOURCED: set by judgement in 4.1.0,
# when the width moved from points to degrees (audit item C5).
BG_MIN_HALF_WIDTH = 0.8
BG_MAX_HALF_WIDTH = 6.0
BG_WIDTH_FACTOR = 2.5            # half-width = factor x 75th-percentile FWHM of the strongest peaks

# UNSOURCED. A relative uncertainty of 3 % of the instrument Gaussian and
# Lorentzian widths, added in quadrature to the uncertainty from the
# calibration scatter (instrument.py, sigma_G and sigma_L). It is added in
# every case, not only when the calibration residuals are smaller.
INST_REL_UNC_FLOOR = 0.03

# UNSOURCED (4.1.1). A standard peak whose η sits at 0 or 1 is censored and
# its fitted error (about 1e-5) means nothing, so it is given this uncertainty
# and so down-weighted in the η trend.
ETA_BOUND_UNC = 0.15        # uncertainty assigned to a standard peak whose eta hit 0 or 1 (censored)
