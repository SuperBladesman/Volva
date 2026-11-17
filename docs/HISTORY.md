# How this history was rebuilt

Völva was not under version control before this repository. Each version was kept as a dated copy on the author's own storage, and the history here was rebuilt from those copies before the first push: one commit for every surviving version, in the order the versions were made, each dated by the last save of its files.

The commits after 4.1.2 prepare the first public release (the licence, the README, the threshold marks, continuous integration, the citation files and this document). They were made before the first push. Two of them change files in `volva/`: one fixes a crash of the command line on Windows when its output goes to a file or a pipe (`volva/cli.py`, with `tests/test_cli.py`), and one adds the sourced and unsourced marks to the comments in `volva/constants.py` without changing any value. `pyproject.toml` now requires Python 3.11. CHANGELOG.md lists them all under 4.1.3.

Nothing in the repository was public before it was first pushed to GitHub, on 16th November 2026. Its open development, issues, pull requests and releases date from then.

## Sources

| Version | Surviving copy | Last saved |
|---|---|---|
| 1.1 | `Volva 1.1 - Crystallite Size and Strain Analysis.ipynb` (description only) and `Völva - Crystallite Size and Strain Analysis Tool.ipynb` (description and code) | 19th August 2025 and 6th November 2025 |
| volva-xrd | the `volva-xrd` folder, a package | 2nd to 5th September 2025 |
| 3.0 | `Volva_3.ipynb` | 6th November 2025 |
| 3.1 | `Volva_3.1.ipynb` | 7th November 2025 |
| 4.0 | `Volva 4.0.ipynb` | 24th September 2026, 09:57 |
| Audit of 4.0 | `Volva_4_0_audit_tests.zip`, the test scripts of the audit | 24th September 2026, 10:09 |
| 4.1.2 | the `Volva_4_1` folder, a package | 24th September 2026, 16:28 |

Apart from 3.1, an edit of 3.0, each version was a rewrite of the one before: the 1.1 desktop app, the volva-xrd package, the 3.0 notebook (which reuses the package's class structure), the three-cell 4.0 notebook and the 4.1 package share little code with their predecessors. The diffs between them are therefore large, and say less than the commit messages.

## How the commits are dated

Each commit carries the time its last file was saved, with one exception.

The 1.1 code survives only in a copy saved on 6th November 2025, three minutes after 3.0. Its description cell is identical to the 1.1 notebook of 19th August 2025, which has that cell and nothing else, and it describes the same desktop app with two methods (3.0 has four). The commit is dated to the 1.1 notebook, and its message says that the code may include changes made up to 6th November.

## The version numbers

The volva-xrd package of September 2025 calls itself 4.0.0 in its `__init__.py`. It came between 1.1 and 3.0, and it is not the 4.0 notebook of September 2026. No copy numbered 2 survives.

## What is missing

- Versions before 1.1 were not kept.
- 4.1.0 and 4.1.1 were made on 24th September 2026, the same day as 4.1.2, and were not saved separately. `CHANGELOG.md` describes all three.
- The written report of the audit of 4.0, whose item numbers (C1–C5, M1–M8) the 4.1.0 changelog entry uses, is not in the repository. Its test scripts are, in `audit/4.0/` in the commit before 4.1.2.

## What was changed for release

Every notebook in the history has been put through the same normalisation, and nothing else:

- outputs and execution counts cleared, because saved outputs contain plots of whatever data a copy was last run on;
- cell metadata and cell ids removed, and notebook metadata reduced to a standard Python 3 kernel (notebook format 4.4);
- one stable file name, `Volva.ipynb`, in place of a new name for each version, so that `git log -p Volva.ipynb` shows what changed between the notebook versions. The original file name is given in each commit message.

From the two package folders, files that were not part of the package were left out: caches (`__pycache__`, `.pytest_cache`), the Google Drive `desktop.ini` files, editor settings (`.idea`), a virtual environment, and a backup copy of a GUI module (`volva-xrd`'s `interfaces/gui_backup.py`). The audit scripts were extracted from their zip file unchanged.

Two files in 4.1.2 were edited, and the commit says so:

- in `CHANGELOG.md`, five passages that gave the composition of four real laboratory patterns, and the date and some numbers of a real LaB₆ calibration, were replaced by general descriptions. The changes they record are unchanged;
- in `tests/test_regression.py`, two sets of parameters that looked rounded from that calibration were replaced by made-up values: the asymmetry in `test_asymmetric_standard_calibrates_smoothly`, now 0.4 exp(−(2θ − 20)/30), and the profile in `test_eta_trend_not_extrapolated_and_bound_peaks_downweighted`, now U = 0.0044, V = −0.0013, W = 0.0039 and η = 0.49 + 0.003 × 2θ. Every test passed before and after the edit.

In rebuilding the history, no code in `volva/`, and nothing in any notebook's code cells, was changed.
