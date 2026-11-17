# Changelog

## Unreleased

### Added
- **Methods text.** Edit ▸ Copy methods text (Ctrl+Shift+C) copies a Methods paragraph describing the analysis as run: the version and its DOI, the peak model and Kα2 handling, the background, the reflections used, the standard and the instrument model, how the instrument was removed, the size convention and K, the six methods, the Monte Carlo and calibration uncertainties, and which method the quoted value comes from, followed by the citation. The *Result* tab shows it, every export writes it as `<name>_methods.txt`, and `volva.report.methods_text(result)` returns it (#32).
- `volva/citation.py` holds the DOI for all versions and the DOI of each release; the About box gives the citation.

### Repository
- CI uses actions/checkout@v5 and actions/setup-python@v6, which run on Node.js 24; GitHub has deprecated Node.js 20 (#22).

## 4.1.3: 16 November 2026

Preparation for the first public release. No change to any result.

### Fixed
- The command line crashed on Windows when its output went to a file or a pipe, which Windows encodes as cp1252: the calibration log and the progress lines contain →, ≈ and η. Output to files and pipes is now UTF-8, and a console shows any character it cannot encode as "?". `tests/test_cli.py` runs the command line under both encodings.

### Repository
- MIT licence (Nik Reeves-McLaren, University of Sheffield).
- Continuous integration: the tests, the offscreen GUI smoke test and the command line on Linux, Windows and macOS, under Python 3.11 to 3.13, and on Windows with the versions pinned in `requirements.txt`.
- `CONTRIBUTING.md`, issue templates, `CITATION.cff` and `.zenodo.json`.
- `docs/HISTORY.md` explains how the history before this release was rebuilt.
- `tests/test_release.py` checks that every copy of the version string agrees; `tools/bump_version.py` sets them for a release.

### Thresholds
- Every default in `constants.py` is marked as sourced (with its reference) or unsourced (with what it was set on), and the README tabulates them. `docs/THRESHOLDS.md` lists, marked the same way, the thresholds still written into the code.
- The comment on `INST_REL_UNC_FLOOR` was wrong: the 3 % is added in quadrature to the calibration uncertainty in every case, not only when the calibration residuals are smaller.

### Packaging
- `requires-python` is now 3.11, the oldest version tested.

### Documentation
- README rewritten for the public repository: installation with pip, a check on the synthetic examples, the command line and the Python interface. Corrected: 4.1.2 had 48 tests, not 42 (there are now 52), and asymmetry is fitted (split pseudo-Voigts, from 4.1.1), so the limitation that said it was not modelled has gone.


## 4.1.2: 24 September 2026

Prompted by slow, noise-dominated analyses of four real laboratory patterns run at a 3σ detection threshold.

### Speed
- Peak fitting now calls `scipy.optimize.least_squares` (trust-region reflective, native bounds) directly. The Kα₂ component is computed inside the model rather than through lmfit constraint expressions, which took most of the run time. Results on the LaB₆ calibration and the synthetic size/strain sweep are identical to 4.1.1 to the printed precision.
- Typical analysis: 15 s → 0.7 s; calibration: 5 s → 0.2 s; the 3σ case that took 1.5–2 min now takes about 1 s. The test suite runs in 20 s rather than 77 s.
- Parallel processing was considered and not adopted: after these changes fitting takes under a second, and starting worker processes in the frozen Windows app would cost more than it saves.

### Peak detection and acceptance
- **Detection significance** is now the net peak area above a local straight baseline, divided by the *local* counting noise (√(gain × expected counts), with the gain estimated from the data). The old test (smoothed prominence against one global σ) accepted hundreds of noise features at 3σ, especially on the high low-angle background and on the flanks of broad peaks. Results are now insensitive to the threshold (3σ and 5σ give the same peaks on the four real patterns). The default is 5σ.
- Features narrower than half the instrument FWHM are ignored at detection. Fitted peaks narrower than 0.7 × the instrument FWHM are rejected as noise or spikes; the fitted width cannot go below that limit.
- Split maxima of one broad peak are merged.
- Within a group, candidates that turn out to be spurious (area < 3σ, width at a bound, narrower than the instrument) are removed and the group is refitted, so noise cannot distort real peaks.
- A peak is accepted if its area is ≥ 5σ and the fit is adequate: R² ≥ threshold, **or** residuals within counting noise (χ²ν ≤ 1.5). This admits weak high-angle reflections on the real patterns that R² alone rejected, which extends the angular range of the linear plots.

### Interpretation
- If the sample η is at 0 or 1 (or undetermined) for most reflections, the Gaussian/Lorentzian split is not determined by the data. Single-line Voigt and Voigt W-H are then shown greyed out and excluded from the strain vote, the anisotropy test and the headline. All four real patterns are in this state: their peaks are more heavy-tailed than a Lorentzian, as a broad or bimodal crystallite-size distribution produces (reproduced with a synthetic bimodal sample).
- With strain present and the shape undetermined, anisotropy is not claimed; the scatter about the Halder-Wagner line is reported as an upper limit.

### Other
- Method table in the text report no longer runs columns together.


## 4.1.1: 24 September 2026

Changes prompted by the first calibration from a real LaB₆ scan.

### Instrument profile
- Peaks are now fitted as **split (asymmetric) pseudo-Voigts**. The real LaB₆ peaks have a strong low-angle tail from axial divergence, largest at low angle and falling with 2θ; symmetric fits turned this into erratic G/L splits. Height, FWHM and integral breadth are independent of the asymmetry parameter.
- **New instrument model:** a Caglioti curve for the *total* FWHM plus a smooth η, with G and L obtained by TCH decomposition of the smooth curves. The 4.1.0 model (separate G and L fits per peak) produced a zig-zag Gaussian width and scattered Lorentzian widths on real data. Old profiles still load and use the old model.
- η is constant unless a linear trend is significant (> 2 SE); peaks whose η hits 0 or 1 are down-weighted; η is not extrapolated outside the calibrated range (previously the instrument Gaussian width could fall to zero above ~105°).
- Fits that stall with η at a bound are retried with the bound parameter fixed.
- New calibration warnings: noisy η, and under-sampled peaks (FWHM < 3 steps). The recommended calibration step is ≤ FWHM/5.

### Interpretation
- Anisotropy is claimed only if it survives dropping the most deviant reflection. If one reflection alone causes the spread (on one of the real patterns, a single weak peak), Völva names it and asks for it to be checked instead of declaring the sample anisotropic.

### Interface
- Instrument tab: total FWHM with the Caglioti curve and derived G/L; η and asymmetry against 2θ; sample/instrument breadth after an analysis.
- Peak inspector, peak CSV and session JSON include the asymmetry (skew).
- The status bar no longer says "Calibrating instrument profile…" after calibration has finished.
- Strain from methods that find no Gaussian broadening reads "none detected".


## 4.1.0: 24 September 2026

This release rebuilds Völva 4.0 as a package with a PySide6 desktop application. The audit IDs (C1–C5, M1–M8) refer to *Völva 4.0 audit for lab rollout*.

### Correctness
- **C1:** pseudo-Voigt FWHM is now 2σ, lmfit's definition. 4.0 overestimated widths by up to 17.7%.
- **C2:** instrumental correction is now Voigt-based: TCH decomposition, with Gaussian widths subtracted in quadrature and Lorentzian widths linearly. The 7–35% size underestimate is removed.
- **C3:** unresolved reflections give "> X nm" bounds, calculated from each reflection's 2σ breadth limit. They are never averaged into a size, in any output.
- **C4:** SSP is removed as a separate method. In its dimensionally consistent form it is identical to Halder-Wagner, and this is now tested. The citation is corrected to Pal (2025).
- **C5:** the SNIP background width is set in degrees and chosen automatically from peak widths, so results no longer depend on scan step size. SNIP is vectorised and works in LLS space.
- **M1:** the Kα₂ component is always fitted when enabled, at its exact Bragg-law position. Resolved Kα₂ lines are recognised and not treated as peaks.
- **M2:** Mo Kα₂ is corrected to 0.0713590 nm. Wavelengths are read from XRDML files.
- **M7:** integral breadths are used throughout. Sizes are L_vol (K = 1), sphere diameter (K = 4/3) or a custom K.
- **M8:** a negative W-H slope is reported as non-physical, not as strain.
- Strain is tested for significance and materiality. Non-significant strain is given as an upper limit.
- Sample Gaussian and Lorentzian components are assigned only if they are significant, so the correction preserves the observed total breadth. Sizes near the resolution limit (340–475 nm synthetic) are now within ±3%.
- If most reflections are near or below the resolution threshold, the headline is a 2σ lower bound and no mean is given.
- The calibration systematic is estimated from ±1σ offset and tan θ tilt of the instrument profile, and strain is re-tested against it.

### Methods
- New methods: single-line Voigt (de Keijser 1982) and Voigt Williamson-Hall (Langford 1992).
- New interpretation layer: the strain vote, an anisotropy χ² test, a Monshi slope test, profile shape, scatter about each line, and method agreement. From these it chooses a headline value and says why.
- Overlapping peaks are fitted together, each with a local linear background.

### Traceability
- **M3:** uncertainties are propagated by Monte Carlo from fit covariances and calibration scatter, with a fixed seed.
- **M4:** exports carry full provenance: the session JSON, SHA-256 hashes of sample and standard files, all settings, library versions, and CSV headers.
- **M5:** a missing instrument profile is flagged as critical and every result is labelled a lower bound.
- **M6:** a mismatch in wavelength or Kα₂ handling between the calibration and the sample is critical, and no values are reported.
- Instrument profiles can be saved and loaded as JSON, with a LaB₆ / Si position check (wavelength and zero error).

### Usability
- The PySide6 application runs from Anaconda Prompt (`python -m volva`) or as a stand-alone PyInstaller build.
- A command-line batch mode (`python -m volva analyse *.xy`) is added.
- A robust loader reads headers, comma decimals, CSV and XRDML.
- A test suite of 42 tests is included, covering end-to-end synthetic recovery and the edge cases found by an independent review.
