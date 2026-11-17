"""The short verdict above the headline says what kind of answer the analysis gives."""
import os

from conftest import sample

from volva.instrument import InstrumentProfile
from volva.io import load_pattern
from volva.pipeline import AnalysisSettings, rerun_methods, run_analysis

EX = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")
LEVELS = {"ok", "info", "warning", "critical"}


def _v(res):
    it = res.interpretation
    assert it.verdict and it.verdict_level in LEVELS
    return it.verdict, it.verdict_level


def test_size_and_strain_resolved():
    r = run_analysis(load_pattern(os.path.join(EX, "SYNTHETIC_CeO2_Lvol34nm_strain0.106pc.xy")),
                     InstrumentProfile.load(os.path.join(EX, "SYNTHETIC_LaB6_profile.json")))
    assert _v(r) == ("Size and strain resolved", "ok")


def test_size_only(profile):
    assert _v(run_analysis(sample(50, seed=11), profile)) == ("Size resolved; no strain resolved", "ok")


def test_no_profile_is_a_lower_bound():
    assert _v(run_analysis(sample(50, seed=11), None)) == ("Lower bound only: no instrument profile", "critical")


def test_beyond_resolution_is_a_lower_bound(profile):
    verdict, level = _v(run_analysis(sample(3000, seed=11), profile))
    assert verdict.startswith("Lower bound only: beyond the instrument's") and level == "warning"


def test_no_peaks_selected(profile):
    r = run_analysis(sample(50, seed=11), profile)
    rerun_methods(r, use=[False] * len(r.peaks))
    assert _v(r) == ("No peaks selected", "warning")


def test_setup_mismatch_is_not_reported(profile):
    r = run_analysis(sample(50, seed=11), profile, AnalysisSettings(doublet=False))
    assert _v(r) == ("Not reported: the settings do not match the instrument profile", "critical")
