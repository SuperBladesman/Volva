# Völva

**Crystallite size and microstrain from powder X-ray diffraction line broadening, with the instrument removed properly and the limits stated.**

[![Licence: MIT](https://img.shields.io/badge/licence-MIT-blue.svg)](LICENSE)
![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)
![Version 4.1.2](https://img.shields.io/badge/version-4.1.2-green.svg)
[![tests](https://github.com/SuperBladesman/Volva/actions/workflows/tests.yml/badge.svg)](https://github.com/SuperBladesman/Volva/actions/workflows/tests.yml)

Most crystallite sizes in the literature come from the Scherrer equation, with the instrument's contribution taken off by subtracting full widths in quadrature. For Voigt-shaped peaks that shortcut is wrong. Völva's own 4.0 notebook used it, and on synthetic patterns it underestimated 20–200 nm crystallites by 7–35 %; it also reported 2 µm crystallites as 124 nm. A size beyond the resolution of the instrument is not a size.

Völva fits every reflection as a pseudo-Voigt Kα₁ + Kα₂ doublet, removes a measured instrument profile component by component (Gaussian widths in quadrature, Lorentzian widths linearly) and runs six methods on what is left. It then says where the methods agree, which of their assumptions the data contradict, and when the honest answer is only a bound. It is a desktop application (PySide6), a command line and a Python package.

It is written for users of laboratory diffractometers who need a size and a strain they can defend, and for the people who check their numbers.

> **Validation status.** Völva recovers known sizes and strains from synthetic patterns to within a few per cent (see [Tests](#tests)). It has not yet been checked against a certified size standard; the protocol, for NIST SRM 1979 and SRM 660c, is in [VALIDATION.md](VALIDATION.md). Until that is done, treat results as provisional.

## What Völva does, and why

**Instrumental correction is Voigt-based.** Each fitted pseudo-Voigt (FWHM, η) is split into Gaussian (G) and Lorentzian (L) widths using the Thompson–Cox–Hastings relations. The instrument profile is modelled as Rietveld programs do: the total FWHM of the standard follows a Caglioti curve, FWHM² = U tan²θ + V tan θ + W, and η is constant (or linear in 2θ if the trend is significant); instrument G and L at any angle come from TCH decomposition of these smooth curves. Standard and sample peaks are fitted as split pseudo-Voigts, so axial-divergence asymmetry does not bias the widths. Sample widths follow: G in quadrature, L linearly.

**Breadths are integral breadths.** With K = 1, the Scherrer equation then gives the volume-weighted column length L_vol, the quantity certified for NIST SRM 1979. *Sphere diameter* (K = 4/3) and a custom K are also available.

**Uncertainties are propagated.** Fit covariances and the calibration scatter are carried through by Monte Carlo, with a fixed seed so that identical inputs give identical outputs. Linear methods report regression standard errors; per-reflection methods report the spread between reflections.

**Calibration error is treated as a systematic.** An error in the instrument profile is common to every reflection, so neither statistic above captures it. Völva re-runs every method with the instrument profile shifted by ±1σ, both uniformly and as a tilt in tan θ (which can mimic strain), and the headline quotes the resulting systematic uncertainty. A strain is called resolved only if it exceeds 2 × (statistical ⊕ calibration) uncertainty.

**Sample components need to be significant.** A Gaussian or Lorentzian sample component is assigned only if the observed width exceeds the instrument's by more than 2σ; otherwise the whole excess breadth goes to the other component. This stops noise in the quadrature subtraction from inflating breadths near the resolution limit.

**Results are honest about limits.**
- A reflection whose sample breadth is below 2σ gives only a lower bound, from its own 2σ breadth limit.
- If most reflections are within the resolution threshold (β_obs/β_inst < 1.2, adjustable), the headline is a 2σ lower bound ("> X nm") and no mean is reported.
- Völva warns when a size exceeds half the largest size the instrument can resolve.
- Strain that is not significant is reported as an upper limit.
- Without an instrument profile, every result is labelled a lower bound.
- If the radiation or Kα₂ setting does not match the instrument profile, no values are reported.

### The six methods and what their contrasts mean

| Method | Size from | Strain from | Assumption |
|---|---|---|---|
| Scherrer | total sample breadth, each reflection | – | all broadening is size |
| Single-line Voigt (de Keijser 1982) | Lorentzian part, each reflection | Gaussian part, each reflection | size is Lorentzian, strain is Gaussian |
| Monshi–Scherrer (2012) | intercept of ln β vs ln(1/cos θ) | slope ≠ 1 signals strain | size only |
| Williamson–Hall (UDM) | intercept | slope | size and strain breadths add linearly |
| Voigt Williamson–Hall (Langford 1992) | Lorentzian intercept | Gaussian slope | size is Lorentzian, strain is Gaussian |
| Halder–Wagner (1966) | slope | intercept | parabolic Voigt approximation |

The *Result* tab turns these contrasts into explicit statements:

- **Strain present?** Völva asks whether the Williamson–Hall, Voigt Williamson–Hall, Halder–Wagner and single-line strains are both significant (above 2σ) and material (at least 15 % of the high-angle breadth), and whether Scherrer sizes fall with sin θ. If strain is present, Scherrer sizes are lower bounds and the headline moves to a size–strain method.
- **Anisotropic?** Per-reflection sizes are tested for spread beyond their uncertainties (χ² test, CV > 20 %). If one reflection alone causes the spread, Völva names it rather than calling the sample anisotropic. If the sizes do differ, a single size is not representative; Völva says so and reports the range. Add hkl labels to see which directions are short.
- **Monshi slope** consistent with 1 means size-only broadening.
- **Profile shape.** A mostly Lorentzian sample profile points to size broadening, a mostly Gaussian one to strain or defects. If η sits at 0 or 1 for most reflections the split is not determined, and the two methods that depend on it are left out of the decision.
- **Agreement.** The spread of sizes across methods is shown against their differing assumptions.

The size–strain plot (SSP) is not reported separately. In its dimensionally consistent form it is algebraically identical to Halder–Wagner (Pal, *J. Alloys Compd.* **1048**, 185324, 2025); the test `test_ssp_dimensionally_consistent_form_equals_hw` checks this.

## Getting started

You need Python 3.11 or later. Völva is tested with Python 3.11 to 3.13 on Linux, Windows and macOS.

```bash
pip install "git+https://github.com/SuperBladesman/Volva.git"
```

This installs two commands: `volva`, the desktop application, and `volva-cli`, the command line. `python -m volva` does the same as `volva-cli`, and starts the application if given no arguments.

### Check the installation

The repository has synthetic example patterns (`examples/`, described in [examples/README.md](examples/README.md)). From a copy of the repository:

```bash
git clone https://github.com/SuperBladesman/Volva.git
cd Volva
volva-cli analyse examples/SYNTHETIC_CeO2_Lvol34nm.xy --profile examples/SYNTHETIC_LaB6_profile.json -o results
```

The headline should be L_vol = 33.9 ± 0.2 nm from Scherrer over 13 reflections, with strain not resolved; the pattern was made with L_vol = 33.9 nm and no strain. To run the tests as well:

```bash
pip install ".[test]"
python -m pytest
```

### Command line

```bash
volva-cli calibrate LaB6_scan.xy -o my_instrument_Cu.json --note "Cu, 1/4° divergence slit" --figure
volva-cli analyse data/*.xy --profile my_instrument_Cu.json -o results
```

`volva-cli calibrate --help` and `volva-cli analyse --help` list the options, which match the settings below.

### From Python

```python
from volva import load_pattern, InstrumentProfile, run_analysis

profile = InstrumentProfile.load("examples/SYNTHETIC_LaB6_profile.json")
result = run_analysis(load_pattern("examples/SYNTHETIC_CeO2_Lvol34nm_strain0.106pc.xy"), profile)

print(result.interpretation.headline)
for key, m in result.methods.items():
    print(m.name, m.size_nm, m.strain)
```

`calibrate_from_standard` builds a profile from a standard pattern, `AnalysisSettings` has the settings, and `rerun_methods` repeats the methods after peaks are deselected or labelled.

### Windows: a conda environment, or a stand-alone app

With Anaconda, from Anaconda Prompt in a copy of the repository:

```bat
conda env create -f environment.yml
conda activate volva41
python -m pytest -q
python -m volva
```

`launch_volva.bat` activates the environment and starts the application, so users can double-click it. If Anaconda is not in `%USERPROFILE%\anaconda3`, edit `CONDA_ROOT` in that file.

For PCs without Anaconda, build a stand-alone app once, on any Windows PC that has it:

```bat
conda activate volva41
build_windows.bat
```

This runs the tests and builds the app with PyInstaller, which takes a few minutes. Output goes to `%USERPROFILE%\Volva_build\`, outside any synced folder, because it is about 350 MB:

- `dist\Volva\Volva.exe` is the application. Double-click it to check it works.
- `dist\Volva\Create desktop shortcut.bat` adds Völva to the desktop and the Start menu.
- `Volva_4.1.2_Setup.exe` is a single-file installer, built only if the free [Inno Setup 6](https://jrsoftware.org/isinfo.php) is installed. It installs per user, needs no administrator rights, and creates desktop and Start-menu shortcuts and an uninstaller.

To put it on another PC, run the installer there, or copy the whole `dist\Volva` folder (*e.g.* to `C:\Users\<name>\AppData\Local\Programs\Volva`) and double-click *Create desktop shortcut.bat* inside it. The executable is not code-signed, so Windows SmartScreen may say "Windows protected your PC": choose *More info → Run anyway*. On managed PCs, IT may need to allow it. The first start takes 10–20 s.

To update, rebuild, then reinstall or replace the folder. Instrument profiles (`.json`) are stored wherever users save them, not inside the app, so they survive updates.

## Workflow

1. **Instrument profile (once per instrument configuration).** Measure NIST SRM 660c LaB₆ with exactly the optics, slits and scan programme used for samples. Use **Calibrate…**, check the **Instrument** tab and save the profile (`.json`). Völva reloads the last profile at start-up. Recalibrate after any change to the optics, or at the lab's agreed interval.
2. **Open** the sample pattern: `.xy`, `.xye`, `.txt`, `.dat`, `.csv`, `.asc` or `.xrdml`. For XRDML files, the wavelength recorded by the instrument is used.
3. **Analyse.** Völva fits the background, detects peaks, fits every peak as a pseudo-Voigt Kα₁ + Kα₂ doublet (overlapping peaks together), removes the instrument profile and runs all six methods.
4. **Check the fits.** Use the *Pattern*, *Peak* (◀ ▶) and *Gallery* tabs. Untick doubtful peaks in the table, and type hkl labels if known; the results update straight away.
5. **Read the *Result* tab.** It gives a headline value, what the methods say together, and the full method table.
6. **File ▸ Export results** writes a report, CSV tables, figures and a session file (see [Outputs](#outputs)).

## Outputs

| File | Contents |
|---|---|
| `<name>_report.txt` | headline, findings, method table, per-reflection table, method notes, instrument profile, full log |
| `<name>_peaks.csv` | per reflection: position, FWHM, η, G/L (observed, instrument, sample), breadths, ratio, Scherrer size or bound, flags |
| `<name>_summary.csv` | per method: size ± error or bound, strain ± error or limit, R², n, warnings |
| `<name>_session.json` | complete record: software and library versions, SHA-256 of sample and standard files, all settings, profile parameters, results, log |
| `<name>_*.png` | pattern, sizes by reflection, comparison, Williamson–Hall, Voigt Williamson–Hall, Halder–Wagner, instrument (300 dpi) |

Every CSV starts with `#` header lines giving the version, run time, file checksum, instrument profile, wavelength and K.

## Settings

| Setting | Default | Notes |
|---|---|---|
| Radiation | Cu Kα | *From file* for XRDML; *Custom* for other sources |
| Fit Kα₂ | on | Turn off only for monochromatic (Kα₁) data. It must match the calibration; Völva refuses to report if not |
| Size reported as | L_vol (K = 1) | Sphere diameter (K = 4/3), or a custom K |
| Background half-width | automatic | SNIP clipping width in **degrees**, set from peak widths. Völva warns if it is smaller than 1.5 × the widest FWHM |
| Detection threshold | 5σ | Net peak area above a local baseline, in standard deviations of the local counting noise. Features narrower than the instrument allows are ignored |
| Minimum fit R² | 0.90 | Fit quality needed to accept a peak. A weak peak with lower R² is still accepted if its residuals are within counting noise (χ²ν ≤ 1.5). Every accepted peak also needs a fitted area of at least 5σ and a width of at least 0.7 × the instrument FWHM |
| Resolution limit | 1.2 | β_obs/β_inst below this flags a reflection as near the resolution limit. If most reflections are, only a lower bound is reported |

## Thresholds and where they come from

Every default in `volva/constants.py` is marked there as sourced, with its reference, or unsourced, with what it was set on. This table matches it. Thresholds still written into the code are listed, marked the same way, in [docs/THRESHOLDS.md](docs/THRESHOLDS.md).

| Constant | Value | What it decides | Source |
|---|---|---|---|
| `DEFAULT_MIN_R2` | 0.90 | Minimum R² to accept a sample peak | Unsourced: 0.8 in 3.x, 0.65 in 4.0, 0.90 from 4.1; no measurement recorded |
| `DEFAULT_STD_MIN_R2` | 0.95 | Minimum R² to accept a standard peak | Unsourced: stricter than for samples because standard peaks are strong (4.0 used 0.8) |
| `DEFAULT_PROMINENCE_SIGMA` | 5σ | Detection: net area above a local baseline over the local counting noise | Unsourced: 4.0's 5σ height test, redefined in 4.1.2; on four real patterns 3σ and 5σ selected the same peaks |
| `MIN_AMPLITUDE_SIGNIFICANCE` | 3 | A fitted peak with area below 3 standard errors is spurious and removed from its group | Unsourced (4.1.2) |
| `MIN_AREA_SIGNIFICANCE_ACCEPT` | 5 | An accepted peak needs a fitted area of at least 5 standard errors | Unsourced: equal to the detection default |
| `MIN_FWHM_FRACTION_OF_INSTRUMENT` | 0.7 | A sample peak narrower than 0.7 × the instrument FWHM is noise or a spike | Unsourced: the physical limit is 1; no measurement recorded for the margin |
| `MAX_CHI2_RED_WEAK` | 1.5 | A fit with χ²ν at or below this is adequate even if R² is low | Unsourced: residuals within counting noise give χ²ν ≈ 1 |
| `DEFAULT_MIN_REL_HEIGHT` | 0.01 | Detection floor, as a fraction of the strongest peak | Unsourced: 4.0's default threshold, kept as a floor |
| `DEFAULT_RESOLUTION_RATIO` | 1.2 | β_obs/β_inst below this: near the resolution limit | Unsourced: introduced in 3.1 (November 2025) as "sample width < 20 % of instrumental" |
| `DEFAULT_MC_SAMPLES` | 2000 | Monte Carlo draws for the uncertainties | Unsourced: no convergence study recorded |
| `MIN_PEAKS_LINEAR`, `MIN_PEAKS_MONSHI` | 3 | Fewest reflections for the linear methods | Unsourced: the fewest that leave one degree of freedom for a straight line |
| `BG_MIN_HALF_WIDTH`, `BG_MAX_HALF_WIDTH`, `BG_WIDTH_FACTOR` | 0.8°, 6.0°, 2.5 | SNIP half-width: 2.5 × the 75th-percentile FWHM, within 0.8–6.0° | SNIP is sourced (Ryan *et al.* 1988); the limits and factor are unsourced (4.1.0) |
| `INST_REL_UNC_FLOOR` | 0.03 | 3 % of the instrument G and L widths, added in quadrature to the calibration uncertainty | Unsourced |
| `ETA_BOUND_UNC` | 0.15 | Uncertainty given to a standard peak's η at 0 or 1, so it is down-weighted in the η trend | Unsourced (4.1.1) |

The reference values in the same file are sourced: X-ray wavelengths from *International Tables for Crystallography* Vol. C, Table 4.2.2.1 (the Kα₂/Kα₁ intensity ratio of 0.5 is the usual approximation, not a value from the table); L_vol with K = 1 from Stokes and Wilson (1942) and the sphere conversion from Langford and Wilson (1978); the LaB₆ lattice parameter from the SRM 660c certificate. The silicon lattice parameter is nominal; it is used only for the position check.

## Tests

`python -m pytest` runs 52 tests, in under a minute on a laptop:

- **Profile mathematics:** FWHM = 2σ for lmfit's pseudo-Voigt; TCH round trip; Voigt integral breadth against numerical integration; split pseudo-Voigt area and FWHM.
- **File formats:** headers, comma decimals, delimiters, XRDML with wavelength.
- **Method algebra:** exact recovery from noise-free breadths; SSP ≡ Halder–Wagner.
- **End-to-end recovery** from synthetic patterns (true Voigt convolution, Kα doublet, Poisson noise).
- **Command line:** calibration and analysis of the examples, with output encoded as UTF-8 and as cp1252 (what Windows gives a pipe or a file).
- **Release:** every copy of the version string agrees.

`tests/gui_smoke.py` drives the application without a display (`QT_QPA_PLATFORM=offscreen`): it loads the example profile and pattern, runs an analysis, deselects a peak and re-runs. Continuous integration runs both on Linux, Windows and macOS.

The synthetic sweep made for 4.1 gave these results (the tests themselves check looser tolerances):

| Case | Result |
|---|---|
| Strain-free, 10–120 nm | Scherrer L_vol within −0.4 % to +2.1 % |
| e = 0.05–0.15 %, 30–120 nm | size within −9 % to 0 %; strain within −10 % to +2 % (Voigt Williamson–Hall) |
| Strain-free, 340–475 nm (near the resolution limit) | within −2 % to +1 %, with a calibration systematic of ± 30–70 nm and a measurable-range warning |
| ≥ about 600 nm on this synthetic instrument | reported only as a 2σ lower bound |
| Step size 0.01° / 0.02° / 0.03° (5 nm crystallites) | results agree within 6 % |
| 1–5 µm crystallites | reported as "> X nm" (always below the true size), never as a number |
| Other checks | anisotropy detected when real, not flagged for noisy data; missing profile flagged; wrong wavelength caught by the LaB₆ position check; Kα₂ or radiation mismatch blocks the result; single-reflection and odd-setting edge cases; identical input gives identical output |

Two caveats apply. The synthetic data use the same Voigt model that Völva assumes, so validation on real standards ([VALIDATION.md](VALIDATION.md)) is still needed. Small strain in very small crystallites (10 nm, e ≤ 0.15 %) is often not resolved; that is the detection limit of integral-breadth methods, and Völva reports it as mixed or unresolved rather than guessing.

## Known limitations

- Not yet validated against a certified size standard (see the status note above).
- Integral-breadth methods cannot separate a size distribution from crystallite shape, or treat stacking faults explicitly. For those, use whole-powder-pattern modelling (*e.g.* WPPM, or TOPAS or GSAS-II with microstructure models). Heavy-tailed, super-Lorentzian peaks from a broad size distribution are detected and reported, but not modelled.
- Strain is isotropic in all the linear methods. Anisotropic strain shows up as scatter and is reported as such, not modelled.
- Asymmetry is fitted (split pseudo-Voigts), but reflections below about 20° 2θ, where axial divergence is strongest, deserve extra scrutiny.
- Binary formats (.raw, .brml) are not read. Export the scan as .xy or .xrdml.

Each is, or will be, an open issue.

## Versions

The current version is **4.1.2**. [CHANGELOG.md](CHANGELOG.md) describes each change from 4.1.0, the rebuild of the 4.0 notebook as a package. The history in this repository before the first public release was rebuilt from dated copies of earlier versions, from the 1.1 desktop app of August 2025; [docs/HISTORY.md](docs/HISTORY.md) explains how, and what is missing.

## Citing Völva

Please cite the version you used. Each tagged release will be archived on Zenodo with its own DOI; until the first one, cite the repository and the version (`volva-cli --version`). [CITATION.cff](CITATION.cff) has the details in machine-readable form, and GitHub's "Cite this repository" button reads it.

## Contributing

Questions, bug reports and patterns Völva gets wrong are welcome as GitHub issues. [CONTRIBUTING.md](CONTRIBUTING.md) says what to include, and how to change the code.

## References

- Caglioti, Paoletti & Ricci (1958) *Nucl. Instrum.* 3, 223.
- de Keijser, Langford, Mittemeijer & Vogels (1982) *J. Appl. Cryst.* 15, 308.
- Halder & Wagner (1966) *Acta Cryst.* 20, 312.
- Langford (1992) *NIST Spec. Publ.* 846, 110.
- Langford & Wilson (1978) *J. Appl. Cryst.* 11, 102.
- Monshi, Foroughi & Monshi (2012) *World J. Nano Sci. Eng.* 2, 154.
- Pal (2025) *J. Alloys Compd.* 1048, 185324.
- Ryan *et al.* (1988) *Nucl. Instrum. Methods B* 34, 396 (SNIP).
- Stokes & Wilson (1942) *Proc. Camb. Phil. Soc.* 38, 313.
- Thompson, Cox & Hastings (1987) *J. Appl. Cryst.* 20, 79.
- Williamson & Hall (1953) *Acta Metall.* 1, 22.
- NIST SRM 660c (LaB₆) and SRM 1979 (ZnO) certificates.

## Licence

MIT, © Nik Reeves-McLaren, University of Sheffield. See [LICENSE](LICENSE).
