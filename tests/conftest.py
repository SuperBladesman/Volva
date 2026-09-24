import os
import sys
import warnings

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
warnings.filterwarnings("ignore")

import synth  # noqa: E402
from volva.constants import WAVELENGTHS  # noqa: E402
from volva.instrument import calibrate_from_standard  # noqa: E402
from volva.io import Pattern  # noqa: E402

CU = WAVELENGTHS["Cu Ka"]
TT = np.arange(10, 125, 0.02)


def lvol_true(D):
    """Synthetic generator defines D via Lorentzian FWHM with K=0.94; convert to integral-breadth L_vol."""
    return D * 1.0 / (0.94 * np.pi / 2)


@pytest.fixture(scope="session")
def profile():
    t, y, _ = synth.pattern(synth.LAB6(), None, 0, tt=TT, seed=101)
    prof, log, fits, bg = calibrate_from_standard(Pattern(t, y, path="lab6.xy"), CU)
    assert prof is not None
    return prof


def sample(D, eps=0.0, seed=1, tt=TT, **kw):
    t, y, _ = synth.pattern(synth.CEO2(), D, eps, tt=tt, seed=seed, **kw)
    return Pattern(t, y, path=f"ceo2_{D}.xy")
