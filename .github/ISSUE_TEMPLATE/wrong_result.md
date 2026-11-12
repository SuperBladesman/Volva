---
name: A result that looks wrong
about: A size, strain, bound or warning that does not look right for a pattern
---

**Version.** `volva-cli --version`, or Help ▸ About Völva in the application; installed with pip, conda or the Windows installer?

**Python and system.** `python -c "import sys, platform, numpy, scipy, lmfit; print(sys.version, platform.platform(), numpy.__version__, scipy.__version__, lmfit.__version__)"`

**The pattern and the profile.** The `_session.json` file from **File ▸ Export results** records the settings, the instrument profile and the checksums of both files, so please attach it. Attach the sample and standard patterns too if you are free to share them; if you are not, say what the sample is and roughly what size you expected.

**What Völva reported.** The headline, and any warnings, from the *Result* tab or `_report.txt`.

**What you expected, and why.** For example, a size from another method or program, or a certified value.
