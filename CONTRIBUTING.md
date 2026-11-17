# Contributing to Völva

Völva is maintained by Nik Reeves-McLaren at the University of Sheffield, in the time that research and teaching leave. Issues and pull requests are read and answered as time allows; there is no guaranteed response time.

## Reporting a result that looks wrong

A pattern that Völva gets wrong is the most useful report there is. Many of the changes in the changelog began as one, and each became a test. Please open an issue with:

- the version (`volva-cli --version`, or Help ▸ About Völva in the application) and how you installed it;
- your Python and operating system;
- the `_session.json` file from **File ▸ Export results**, which records the settings, the instrument profile and checksums of the files;
- the patterns themselves if you are free to share them. If you are not, say what the sample is, and what size you expected and why;
- what Völva reported (the headline and any warnings) and what you expected.

Questions about what a result means are welcome as issues too. The README's [What Völva does, and why](README.md#what-völva-does-and-why) and the *Result* tab answer most of them.

## Asking for help

Open an issue. If you are not sure whether something is a bug, open one anyway and say so.

## Changing the code

1. Open an issue first for anything beyond a small fix, so that the change can be agreed before you write it.
2. Fork the repository and work on a branch. Install with `pip install -e ".[test]"`.
3. Add a test that fails before your change and passes after it. The tests use synthetic patterns made by `tests/synth.py`, so that the answer is known; add a case there rather than a measured pattern.
4. Run `python -m pytest`, and `QT_QPA_PLATFORM=offscreen python tests/gui_smoke.py` if you changed the application.
5. A new threshold or default goes in `volva/constants.py`, marked as sourced (with its reference) or unsourced (with what it was set on), and in the README table. A threshold you change keeps its mark up to date.
6. Describe the change in `CHANGELOG.md` under *Unreleased*: what was wrong, what changed, and what it does to results.
7. Open a pull request that names the issue. Continuous integration runs the tests on Linux, Windows and macOS.

Never commit data that cannot be redistributed: measured patterns from other people's samples, or files from commercial databases. Synthetic patterns are fine.

## Releases

Releases are made once a month while Völva is in active development. Each is tagged (`v4.1.3`) and published as a GitHub release, which Zenodo archives with its own DOI. Before tagging, run `python tools/bump_version.py X.Y.Z YYYY-MM-DD`, which sets the version in `pyproject.toml`, `volva/__init__.py`, `CITATION.cff` (with `date-released`), `.zenodo.json`, the Windows build files and the README, and dates the changelog entry; `tests/test_release.py` checks that they agree. After Zenodo has archived the release, add its DOI to `CITATION.cff`, the README and `VERSION_DOIS` in `volva/citation.py`, so that the methods text cites that version.

## Licence

By contributing you agree that your contribution is released under the [MIT licence](LICENSE).
