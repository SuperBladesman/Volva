import csv
import json

from conftest import sample

from volva import __version__
from volva.instrument import InstrumentProfile
from volva.pipeline import run_analysis
from volva.report import write_bundle


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
