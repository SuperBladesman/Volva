import csv
import json

from conftest import sample

from volva import __version__
from volva.instrument import InstrumentProfile
from volva.pipeline import run_analysis
from volva.citation import CONCEPT_DOI, VERSION_DOIS, citation, doi_for
from volva.pipeline import AnalysisSettings
from volva.report import methods_text, write_bundle


def test_bundle_and_provenance(tmp_path, profile):
    r = run_analysis(sample(50, seed=11), profile)
    paths = write_bundle(r, str(tmp_path), figures=True, dpi=80)
    names = {p.split("/")[-1].split("\\")[-1] for p in paths}
    assert any(n.endswith("_session.json") for n in names)
    assert sum(n.endswith(".png") for n in names) == 7
    js = json.load(open([p for p in paths if p.endswith("_session.json")][0], encoding="utf-8"))
    assert js["provenance"]["software"] == f"Völva {__version__}"
    assert js["provenance"]["instrument_profile"]["U"] == profile.U
    assert js["provenance"]["settings"]["seed"] == r.settings.seed
    csvp = [p for p in paths if p.endswith("_summary.csv")][0]
    head = open(csvp, encoding="utf-8").read().splitlines()[0]
    assert __version__ in head
    rows = [row for row in csv.reader(line for line in open(csvp, encoding="utf-8") if not line.startswith("#"))]
    assert rows[0][0] == "method" and len(rows) == 7


def test_profile_roundtrip(tmp_path, profile):
    p = tmp_path / "prof.json"
    profile.save(str(p))
    q = InstrumentProfile.load(str(p))
    assert q.U == profile.U and q.X == profile.X and q.peaks == profile.peaks
    assert abs(q.G_at(60.0) - profile.G_at(60.0)) < 1e-12


def test_methods_text_follows_the_settings(profile):
    pat = sample(50, seed=11)
    base = methods_text(run_analysis(pat, profile))
    assert f"Völva {__version__}" in base and "doi:10.5281/zenodo." in base
    assert "Kα1/Kα2 doublet" in base and "K = 1" in base and "LaB6 (NIST SRM 660c)" in base
    assert "2000 draws" in base and "±1σ" in base and "https://doi.org/" in base
    assert "—" not in base                      # no em dashes in text meant for a manuscript
    sphere = methods_text(run_analysis(pat, profile, AnalysisSettings(size_convention="Sphere diameter (K = 4/3)")))
    assert "sphere diameters" in sphere and "K = 4/3" in sphere
    custom = methods_text(run_analysis(pat, profile, AnalysisSettings(k_custom=0.9)))
    assert "K = 0.9" in custom
    mc = methods_text(run_analysis(pat, profile, AnalysisSettings(mc_samples=500)))
    assert "500 draws" in mc
    lim = methods_text(run_analysis(pat, profile, AnalysisSettings(two_theta_max=90.0)))
    assert "and 90.0° 2θ were analysed" in lim
    none = methods_text(run_analysis(pat, None))
    assert "every size is a lower bound" in none and "±1σ" not in none and "instrument calibration" not in none
    bare = methods_text(run_analysis(pat, profile), with_citation=False)
    assert "https://doi.org/" not in bare


def test_methods_text_counts_the_reflections_used(profile):
    from volva.pipeline import rerun_methods
    r = run_analysis(sample(50, seed=11), profile)
    n = len(r.peaks)
    assert f"All {n} fitted reflections" in methods_text(r)
    rerun_methods(r, use=[i != 0 for i in range(n)])
    assert f"Of the {n} fitted reflections, {n - 1}" in methods_text(r)


def test_citation_uses_the_version_doi_when_there_is_one():
    for v, (doi, year) in VERSION_DOIS.items():
        assert doi.startswith("10.5281/zenodo.") and doi != CONCEPT_DOI and 2025 <= year
        assert doi_for(v) == (doi, year, True)
        assert f"(version {v})" in citation(v) and doi in citation(v)
    d, _, own = doi_for("0.0.0")
    assert d == CONCEPT_DOI and not own and "no DOI of its own" in citation("0.0.0")


def test_bundle_includes_methods_text(tmp_path, profile):
    r = run_analysis(sample(50, seed=11), profile)
    paths = write_bundle(r, str(tmp_path), figures=False)
    mt = [p for p in paths if p.endswith("_methods.txt")]
    assert len(mt) == 1
    txt = open(mt[0], encoding="utf-8").read()
    assert txt.startswith("Crystallite size and microstrain were determined") and "Zenodo" in txt
