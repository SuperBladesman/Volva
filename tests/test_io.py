import numpy as np
import pytest

from volva.io import LoadError, load_pattern

TT = np.arange(10, 80, 0.02)
Y = 100 + 1000 * np.exp(-((TT - 30) / 0.1) ** 2)


def _w(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return str(p)


@pytest.mark.parametrize("name,maker", [
    ("plain_header.xy", lambda: "2theta intensity\n" + "\n".join(f"{a:.4f} {b:.1f}" for a, b in zip(TT, Y))),
    ("hash.xy", lambda: "# comment\n" + "\n".join(f"{a:.4f} {b:.1f}" for a, b in zip(TT, Y))),
    ("csv_header.csv", lambda: "2Theta,Counts\n" + "\n".join(f"{a:.4f},{b:.1f}" for a, b in zip(TT, Y))),
    ("semicolon_comma_decimal.txt", lambda: "\n".join(f"{a:.4f};{b:.1f}".replace(".", ",") for a, b in zip(TT, Y))),
    ("space_comma_decimal.txt", lambda: "\n".join(f"{a:.4f} {b:.1f}".replace(".", ",") for a, b in zip(TT, Y))),
    ("tab.dat", lambda: "\n".join(f"{a:.4f}\t{b:.1f}" for a, b in zip(TT, Y))),
    ("three_col.xye", lambda: "\n".join(f"{a:.4f} {b:.1f} {np.sqrt(b):.2f}" for a, b in zip(TT, Y))),
    ("trailing.xy", lambda: "\n".join(f"{a:.4f} {b:.1f}" for a, b in zip(TT, Y)) + "\nEND\n"),
    ("bruker_meta.txt", lambda: "[Measurement]\nWavelength=1.5406\n[Data]\n" + "\n".join(f"{a:.4f} {b:.1f}" for a, b in zip(TT, Y))),
])
def test_text_formats(tmp_path, name, maker):
    p = load_pattern(_w(tmp_path, name, maker()))
    assert len(p.two_theta) == len(TT)
    assert abs(p.two_theta[0] - 10) < 1e-6 and abs(p.intensity.max() - Y.max()) < 1
    assert len(p.sha256) == 64


def test_xrdml_with_wavelength(tmp_path):
    counts = " ".join(f"{v:.0f}" for v in Y)
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<xrdMeasurements xmlns="http://www.xrdml.com/XRDMeasurement/1.5">
 <xrdMeasurement>
  <usedWavelength intended="K-Alpha 1">
   <kAlpha1 unit="Angstrom">1.5405980</kAlpha1><kAlpha2 unit="Angstrom">1.5444260</kAlpha2>
   <ratioKAlpha2KAlpha1>0.5</ratioKAlpha2KAlpha1>
  </usedWavelength>
  <scan><dataPoints>
   <positions axis="2Theta" unit="deg"><startPosition>{TT[0]}</startPosition><endPosition>{TT[-1]}</endPosition></positions>
   <commonCountingTime unit="seconds">1.0</commonCountingTime>
   <counts unit="counts">{counts}</counts>
  </dataPoints></scan>
 </xrdMeasurement>
</xrdMeasurements>"""
    p = load_pattern(_w(tmp_path, "scan.xrdml", xml))
    assert len(p.two_theta) == len(TT)
    assert abs(p.wavelength["ka1"] - 0.154059800) < 1e-9
    assert abs(p.wavelength["ka2"] - 0.154442600) < 1e-9


def test_rejects_binary_and_garbage(tmp_path):
    with pytest.raises(LoadError):
        load_pattern(_w(tmp_path, "x.raw", "binary"))
    with pytest.raises(LoadError):
        load_pattern(_w(tmp_path, "x.txt", "hello\nworld\n"))
