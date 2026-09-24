"""End-to-end recovery of known size and strain from synthetic patterns (true Voigt, Kα doublet, Poisson noise)."""
import numpy as np
import pytest
from conftest import CU, TT, lvol_true, sample

from volva.pipeline import AnalysisSettings, run_analysis
from volva.io import Pattern
from volva.instrument import calibrate_from_standard
import synth


@pytest.mark.parametrize("D", [10, 50, 100])
def test_size_recovery_strain_free(profile, D):
    r = run_analysis(sample(D, seed=D), profile)
    truth = lvol_true(D)
    assert abs(r.methods["scherrer"].size_nm - truth) / truth < 0.05
    assert abs(r.methods["hw"].size_nm - truth) / truth < 0.06
    assert r.interpretation.strain_state == "not_resolved"
    assert "Scherrer" in r.interpretation.headline_method


def test_size_and_strain_recovery(profile):
    r = run_analysis(sample(50, 0.001, seed=3), profile)
    truth = lvol_true(50)
    e_true = 0.001 * 0.5 * np.sqrt(np.pi / np.log(2))      # integral-breadth convention
    for key in ("voigt_wh", "hw"):
        m = r.methods[key]
        assert abs(m.size_nm - truth) / truth < 0.08, key
        assert abs(m.strain - e_true) / e_true < 0.2, key
    assert r.interpretation.strain_state == "resolved"
    assert r.methods["scherrer"].size_nm < 0.85 * truth          # Scherrer is a lower bound here
    assert r.methods["monshi"].extra["slope_t"] > 2               # Monshi slope flags strain


@pytest.mark.parametrize("D", [1500, 5000])
def test_resolution_limited_reported_as_bound(profile, D):
    """Völva 4.0 reported 2 µm crystallites as '124 nm'. Now: a lower bound, never a mean."""
    r = run_analysis(sample(D, seed=4), profile)
    h = r.interpretation.headline
    assert ">" in h and "lower bound" in h
    lb = float(h.split(">")[1].split("nm")[0])
    assert lb < lvol_true(D)                                   # a true lower bound
    assert r.interpretation.headline_method == "resolution bound"


def test_step_size_invariance_broad_peaks():
    """
    Völva 4.0 SNIP worked in points; 4.1 works in degrees, so the same reflection must give
    the same size whatever the step. At L_vol ≈ 3.4 nm several CeO2 pairs (311/222, 331/420,
    531/600) overlap so strongly that they appear as one peak without phase information; the
    comparison is therefore made reflection by reflection, and accuracy is checked on the
    isolated reflections.
    """
    from volva.instrument import calibrate_from_standard
    per_step = []
    for step in (0.01, 0.02, 0.03):
        tt = np.arange(10, 125, step)
        t, y, _ = synth.pattern(synth.LAB6(), None, 0, tt=tt, seed=21)
        prof, *_ = calibrate_from_standard(Pattern(t, y, path="s.xy"), CU)
        r = run_analysis(sample(5, seed=22, tt=tt), prof)
        rows = r.methods["scherrer"].per_peak
        per_step.append({round(row["two_theta"]): (row["size"], row["size_err"]) for row in rows
                         if row.get("size") and row.get("size_err")})
    merged = {56, 57, 58, 59, 76, 77, 78, 79, 114, 115, 116, 117}   # unresolvable pairs at this size
    common = set.intersection(*(set(d) for d in per_step)) - merged
    assert len(common) >= 5
    for k in common:
        # each step is a different noise realisation: agreement within 3 combined standard errors
        v = np.array([d[k][0] for d in per_step])
        e = np.array([d[k][1] for d in per_step])
        assert np.ptp(v) < 3 * np.sqrt(2) * np.max(e) and np.ptp(v) / np.mean(v) < 0.15, (k, v, e)
    isolated = [k for k in (29, 33, 47, 88, 95) if k in common]
    assert len(isolated) >= 4
    means = [np.mean([d[k][0] for k in isolated]) for d in per_step]
    assert np.ptp(means) / np.mean(means) < 0.06
    mean = np.mean(means)
    assert abs(mean - lvol_true(5)) / lvol_true(5) < 0.08


def test_no_profile_flagged_as_lower_bound():
    r = run_analysis(sample(50, seed=5), None)
    assert r.interpretation.headline.startswith("LOWER BOUND")
    assert any(f.level == "critical" for f in r.interpretation.findings)


def test_anisotropy_detected(profile):
    refl = synth.CEO2()
    big = {round(p, 3) for i, (p, m) in enumerate(refl) if i % 2 == 0}
    t, y, _ = synth.pattern(refl, lambda pos: 80 if round(pos, 3) in big else 30, 0, tt=TT, seed=6)
    r = run_analysis(Pattern(t, y, path="aniso.xy"), profile)
    assert r.interpretation.anisotropic


def test_noisy_data_not_called_anisotropic(profile):
    t, y, _ = synth.pattern(synth.CEO2(), 50, 0, tt=TT, seed=7, scale=800, bg=30)
    r = run_analysis(Pattern(t, y, path="weak.xy"), profile)
    assert not r.interpretation.anisotropic
    assert abs(r.methods["scherrer"].size_nm - lvol_true(50)) / lvol_true(50) < 0.12


def test_reproducible(profile):
    a = run_analysis(sample(50, seed=8), profile)
    b = run_analysis(sample(50, seed=8), profile)
    assert a.methods["scherrer"].size_nm == b.methods["scherrer"].size_nm
    assert a.methods["scherrer"].size_err == b.methods["scherrer"].size_err


def test_calibration_quality(profile):
    assert profile.profile_model == "total" and profile.f_model == "caglioti" and profile.r2_F > 0.9
    assert profile.position_check["ok"]
    x = np.array([30, 60, 90, 110])
    G, L = synth.inst_profile(x)
    assert np.all(np.abs(profile.G_at(x) - G) / G < 0.08)
    assert np.all(np.abs(profile.L_at(x) - L) / L < 0.15)


def test_wrong_wavelength_flagged_by_position_check():
    t, y, _ = synth.pattern(synth.LAB6(), None, 0, tt=TT, seed=31)
    from volva.instrument import calibrate_from_standard
    prof, *_ = calibrate_from_standard(Pattern(t, y, path="s.xy"), dict(ka1=0.1788965, ka2=0.1792850, ratio=0.5))
    assert prof is None or prof.position_check.get("ok") is False


def test_doublet_mismatch_is_critical(profile):
    r = run_analysis(sample(50, seed=9), profile, AnalysisSettings(doublet=False))
    assert any("Kα2 handling differs" in w for w in r.critical)
    assert r.interpretation.headline.startswith("UNRELIABLE")


def test_custom_wavelength_without_ka2_is_caught(profile):
    """Reviewer case: Custom λ with Kα2 = Kα1 on Cu data gave a spurious strain silently."""
    s = AnalysisSettings(wavelength="Custom", custom_ka1=0.1540562, custom_ka2=0.1540562, custom_ratio=0.0)
    r = run_analysis(sample(30, seed=10), profile, s)
    assert r.interpretation.headline.startswith("UNRELIABLE")


def test_near_limit_size_not_biased_low(profile):
    """Reviewer case: ~340-470 nm sizes were reported 15-30% low with no warning."""
    r = run_analysis(sample(500, seed=12), profile)
    it = r.interpretation
    truth = lvol_true(500)
    if it.headline_method == "resolution bound":
        lb = float(it.headline.split(">")[1].split("nm")[0])
        assert lb < truth
    else:
        s = r.methods["scherrer"].size_nm
        assert abs(s - truth) / truth < 0.15 or any(f.topic == "Measurable range" for f in it.findings)


def test_constrained_correction_preserves_total_breadth():
    from volva.correction import constrained_correction
    from volva.profiles import voigt_integral_breadth as vb
    Gi, Li = 0.05, 0.03
    Go, Lo = np.array([0.04, 0.052]), np.array([0.05, 0.034])     # G deficit / both slightly above
    Gs, Ls = constrained_correction(Go, Lo, Gi, Li)
    assert Gs[0] == 0 and Ls[0] > 0
    assert abs(vb(np.sqrt(Gi**2 + Gs**2), Li + Ls)[0] - vb(Go, Lo)[0]) < 1e-9


def test_single_peak_and_odd_settings_do_not_crash(profile):
    r = run_analysis(sample(50, seed=13), profile, AnalysisSettings(resolution_ratio=1.0))
    from volva.pipeline import rerun_methods
    use = [i == 3 for i in range(len(r.peaks))]
    rerun_methods(r, use=use)
    assert "nan" not in r.interpretation.headline
    assert "single reflection" in r.interpretation.headline


def test_asymmetric_standard_calibrates_smoothly():
    """Axial-divergence asymmetry (strong at low angle) must not distort the instrument widths."""
    asym = lambda tt: 0.4 * np.exp(-(tt - 20) / 30)
    t, y, truth = synth.pattern(synth.LAB6(), None, 0, tt=TT, seed=202, asym=asym)
    prof, log, fits, bg = calibrate_from_standard(Pattern(t, y, path="lab6_asym.xy"), CU)
    assert prof is not None and prof.r2_F > 0.9
    sk = [q["skew"] for q in prof.peaks if q["two_theta"] < 40]
    assert np.mean(sk) > 0.2                       # low-angle tail detected
    x = np.array([30, 60, 90, 110])
    G, L = synth.inst_profile(x)
    F = synth.voigt_fwhm(G, L)
    assert np.all(np.abs(prof.fwhm_at(x) - F) / F < 0.06)
    inside = np.linspace(*prof.two_theta_range, 50)
    assert np.all(prof.G_at(inside) > 0.2 * synth.inst_profile(inside)[0])


def test_eta_trend_not_extrapolated_and_bound_peaks_downweighted():
    from volva.instrument import InstrumentProfile
    p = InstrumentProfile(profile_model="total", FU=0.0044, FV=-0.0013, FW=0.0039, eta_coef=[0.49, 0.003],
                          two_theta_range=(20.0, 100.0))
    assert p.eta_at(140.0) == pytest.approx(p.eta_at(100.0))
    assert p.eta_at(5.0) == pytest.approx(p.eta_at(20.0))
    assert p.G_at(140.0) > 0


def test_split_pvoigt_keeps_area_and_fwhm():
    from volva.peaks import split_pvoigt
    x = np.linspace(-3, 3, 60001)
    for sk in (0.0, 0.4, -0.3):
        y = split_pvoigt(x, amplitude=2.0, center=0.0, sigma=0.1, fraction=0.6, skew=sk)
        g = split_pvoigt(x, amplitude=2.0, center=0.0, sigma=0.1, fraction=0.0, skew=sk)
        assert np.sum(g) * (x[1] - x[0]) == pytest.approx(2.0, rel=1e-4)   # Gaussian: no tail truncation
        half = y.max() / 2
        above = x[y >= half]
        assert above.max() - above.min() == pytest.approx(0.2, abs=2e-3)


def test_single_outlier_reflection_does_not_make_sample_anisotropic(profile):
    refl = synth.CEO2()
    odd = round(refl[3][0], 3)
    t, y, _ = synth.pattern(refl, lambda pos: 8 if round(pos, 3) == odd else 50, 0, tt=TT, seed=33)
    r = run_analysis(Pattern(t, y, path="one_odd.xy"), profile)
    it = r.interpretation
    assert not it.anisotropic
    assert any(f.topic == "Anisotropy" and "only because of" in f.text for f in it.findings)


def test_low_detection_threshold_does_not_fit_noise(profile):
    """At a 3σ threshold on weak, noisy data only real reflections may be detected, and fast."""
    import time
    refl = synth.CEO2()
    t, y, _ = synth.pattern(refl, 40, 0, tt=TT, seed=61, scale=1500, bg=80)
    t0 = time.perf_counter()
    r = run_analysis(Pattern(t, y, path="weak.xy"), profile, AnalysisSettings(prominence_sigma=3))
    assert time.perf_counter() - t0 < 20
    true = np.array([p for p, _ in refl])
    for c in r.candidates:
        assert np.min(np.abs(true - c.pos)) < 0.5, f"spurious candidate at {c.pos:.2f}°"
    for p in r.peaks:                       # nothing narrower than the instrument survives
        assert "narrower_than_instrument" not in p.fit.flags


def test_heavy_tailed_profiles_disable_component_methods(profile):
    """A broad size distribution gives super-Lorentzian peaks (η → 1): the G/L split is then
    undetermined and single-line / Voigt W-H must not drive the interpretation."""
    refl = synth.CEO2()
    t, y1, _ = synth.pattern(refl, 12, 0, tt=TT, seed=71, scale=1.2e4, bg=75)
    _, y2, _ = synth.pattern(refl, 120, 0, tt=TT, seed=72, scale=0.8e4, bg=75)
    r = run_analysis(Pattern(t, y1 + y2, path="bimodal.xy"), profile)
    it = r.interpretation
    assert not it.shape_determined
    assert r.methods["single_line"].extra.get("shape_undetermined")
    assert r.methods["voigt_wh"].extra.get("shape_undetermined")
    assert it.headline_key not in ("single_line", "voigt_wh")
