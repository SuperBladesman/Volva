"""Every place the version is written down agrees."""
import json
import re
from pathlib import Path

import volva

ROOT = Path(__file__).resolve().parents[1]


def _read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def test_versions_agree():
    v = volva.__version__
    assert re.search(rf'^version = "{re.escape(v)}"$', _read("pyproject.toml"), re.M)
    assert re.search(rf"^version: {re.escape(v)}$", _read("CITATION.cff"), re.M)
    assert json.loads(_read(".zenodo.json"))["version"] == v
    assert f'#define AppVer "{v}"' in _read("windows/Volva_installer.iss")
    assert f"badge/version-{v}-green.svg" in _read("README.md")
    assert f"The current version is **{v}**" in _read("README.md")


def test_changelog_has_an_entry_for_this_version_or_unreleased():
    head = [line for line in _read("CHANGELOG.md").splitlines() if line.startswith("## ")][0]
    assert head == "## Unreleased" or head.startswith(f"## {volva.__version__}:")
