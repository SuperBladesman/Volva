"""The command line runs end to end, including where the output encoding is not UTF-8."""
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples"


@pytest.mark.parametrize("encoding", ["utf-8", "cp1252"])
def test_calibrate_and_analyse(tmp_path, encoding):
    """cp1252 is what Windows gives a pipe or a file; →, ≈ and η must not crash it."""
    env = dict(os.environ, PYTHONIOENCODING=encoding, PYTHONPATH=str(ROOT), MPLBACKEND="Agg")
    run = lambda *a: subprocess.run([sys.executable, "-m", "volva", *a], cwd=tmp_path, env=env,
                                    capture_output=True, timeout=300)
    cal = run("calibrate", str(EX / "SYNTHETIC_LaB6_standard.xy"), "-o", "profile.json")
    assert cal.returncode == 0, cal.stderr.decode("utf-8", "replace")
    ana = run("analyse", str(EX / "SYNTHETIC_CeO2_Lvol34nm.xy"), "--profile", "profile.json", "-o", "out")
    assert ana.returncode == 0, ana.stderr.decode("utf-8", "replace")
    assert (tmp_path / "out" / "SYNTHETIC_CeO2_Lvol34nm_report.txt").stat().st_size > 0
    assert "L_vol = 33.9" in ana.stdout.decode("utf-8")
