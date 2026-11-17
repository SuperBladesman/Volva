"""
How to cite Völva, and the Zenodo DOI of each release.

The DOI for all versions resolves to the latest release. Each release's own DOI
is added to VERSION_DOIS in the pull request that follows the release, together
with README.md and CITATION.cff.
"""

from __future__ import annotations

from typing import Tuple

from . import __version__

TITLE = "Völva: crystallite size and microstrain from powder X-ray diffraction line broadening"
AUTHOR = "Reeves-McLaren, N."
CONCEPT_DOI = "10.5281/zenodo.24293731"   # all versions
FIRST_RELEASE_YEAR = 2026

# version -> (Zenodo DOI of that release, year of release)
VERSION_DOIS = {
    "4.1.3": ("10.5281/zenodo.24293732", 2026),
}


def doi_for(version: str = __version__) -> Tuple[str, int, bool]:
    """Return (DOI, year, True) for a released version, or the DOI for all versions with False."""
    if version in VERSION_DOIS:
        doi, year = VERSION_DOIS[version]
        return doi, year, True
    return CONCEPT_DOI, FIRST_RELEASE_YEAR, False


def citation(version: str = __version__) -> str:
    """A reference-list entry for the version used."""
    doi, year, own = doi_for(version)
    note = "" if own else " (DOI for all versions; this version has no DOI of its own)"
    return (f"{AUTHOR} ({year}). {TITLE} (version {version}) [Computer software]. Zenodo. "
            f"https://doi.org/{doi}{note}")


def in_text(version: str = __version__) -> str:
    """How to name the software in running text."""
    doi, year, _ = doi_for(version)
    return f"Völva {version} (Reeves-McLaren, {year}; doi:{doi})"
