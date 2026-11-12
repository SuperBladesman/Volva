"""Set the version everywhere it is written down, for a release.

    python tools/bump_version.py 4.1.3 2026-12-14

Updates pyproject.toml, volva/__init__.py, CITATION.cff (version and
date-released), .zenodo.json, the Windows installer script and build script,
and the README; then dates the "Unreleased" section of CHANGELOG.md.
tests/test_release.py checks that they agree.
"""
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sub(path, pattern, repl, count=1):
    p = ROOT / path
    with open(p, encoding="utf-8", newline="") as f:      # keep CRLF in the Windows files
        s = f.read()
    new, n = re.subn(pattern, repl, s, count=count, flags=re.M)
    if n == 0:
        sys.exit(f"{path}: pattern not found: {pattern}")
    with open(p, "w", encoding="utf-8", newline="") as f:
        f.write(new)


def main(version, date):
    if not re.fullmatch(r"\d+\.\d+\.\d+", version) or not re.fullmatch(r"\d{4}-\d{2}-\d{2}", date):
        sys.exit("usage: bump_version.py X.Y.Z YYYY-MM-DD")
    sub("pyproject.toml", r'^version = "[^"]+"', f'version = "{version}"')
    sub("volva/__init__.py", r'^__version__ = "[^"]+"', f'__version__ = "{version}"')
    sub("CITATION.cff", r"^version: [^\r\n]+", f"version: {version}")
    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    if re.search(r"^date-released:", cff, re.M):
        sub("CITATION.cff", r"^date-released: [^\r\n]+", f"date-released: {date}")
    else:
        sub("CITATION.cff", r"^(version: [^\r\n]+)", rf"\1\ndate-released: {date}")
    z = json.loads((ROOT / ".zenodo.json").read_text(encoding="utf-8"))
    z["version"] = version
    (ROOT / ".zenodo.json").write_text(json.dumps(z, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    sub("windows/Volva_installer.iss", r'^#define AppVer "[^"]+"', f'#define AppVer "{version}"')
    sub("windows/Volva_installer.iss", r"Volva_[\d.]+_Setup\.exe", f"Volva_{version}_Setup.exe")
    sub("build_windows.bat", r"Volva_[\d.]+_Setup\.exe", f"Volva_{version}_Setup.exe", count=0)
    sub("README.md", r"badge/version-[\d.]+-green\.svg", f"badge/version-{version}-green.svg")
    sub("README.md", r"!\[Version [\d.]+\]", f"![Version {version}]")
    sub("README.md", r"Volva_[\d.]+_Setup\.exe", f"Volva_{version}_Setup.exe", count=0)
    sub("README.md", r"The current version is \*\*[\d.]+\*\*", f"The current version is **{version}**")
    d = __import__("datetime").date.fromisoformat(date)
    sub("CHANGELOG.md", r"^## Unreleased$", f"## {version}: {d.day} {d.strftime('%B %Y')}")
    print(f"Version {version}, {date}. Check `git diff`, then run the tests.")


if __name__ == "__main__":
    main(*sys.argv[1:3])
