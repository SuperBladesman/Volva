# Völva 4.1

**Völva 4.1 measures crystallite size and microstrain from powder XRD line broadening. It runs six complementary methods and reports where they agree and where they diverge.** It is a desktop application (PySide6) with a command-line batch mode. It replaces the Völva 4.0 notebook.

> **Validation status.** Völva 4.1 recovers known sizes and strains from synthetic patterns to within a few per cent (see *Tests*). It passed an independent adversarial code review; the issues found were fixed and re-checked. It has **not yet been validated against certified reference materials** on the lab's instruments. Until the protocol in `VALIDATION.md` has been completed, treat results as provisional.

---

## Installation (Windows, Anaconda Prompt)

```bat
cd path\to\Volva_4_1
conda env create -f environment.yml
conda activate volva41
python -m pytest -q          :: optional: 42 tests, about one minute
python -m volva              :: starts the application
```

`launch_volva.bat` activates the environment and starts the application; users can double-click it. If Anaconda is not in `%USERPROFILE%\anaconda3`, edit `CONDA_ROOT` in that file.

### Stand-alone Windows app (no Anaconda needed on the PCs that run it)

Build once, on any Windows PC that has Anaconda (your own is fine):

```bat
conda activate volva41
cd path\to\Volva_4_1
build_windows.bat
```

This runs the tests and then builds the app with PyInstaller, which takes a few minutes. Output goes to `%USERPROFILE%\Volva_build\`, outside Google Drive, because the folder is about 350 MB:

- `dist\Volva\Volva.exe` is the application. Double-click it to check it works.
- `dist\Volva\Create desktop shortcut.bat` adds Völva to the desktop and the Start menu.
- `Volva_4.1.2_Setup.exe` is a single-file installer. It is built only if the free [Inno Setup 6](https://jrsoftware.org/isinfo.php) is installed. It installs per user, needs no admin rights, and creates desktop and Start-menu shortcuts plus an uninstaller.

**To put it on another PC**, either run `Volva_4.1.2_Setup.exe` there, or copy the whole `dist\Volva` folder (e.g. to `C:\Users\<name>\AppData\Local\Programs\Volva`) and double-click *Create desktop shortcut.bat* inside it.

**First launch.** The executable is not code-signed, so Windows SmartScreen may say "Windows protected your PC". Choose *More info → Run anyway*. On managed university PCs, IT may need to allow it. The first start takes 10–20 s; later starts are faster.

**Updating.** Rebuild, then reinstall or replace the folder. Instrument profiles (`.json`) are stored wherever users save them, not inside the app, so they survive updates.

---

## Workflow

1. **Instrument profile (once per instrument configuration).** Measure NIST SRM 660c LaB₆ with exactly the optics, slits and scan programme used for samples. Then use **Calibrate…**, check the **Instrument** tab and save the profile (`.json`). Völva reloads the last profile automatically at start-up. Recalibrate after any change to the optics, or at the lab's agreed interval.
2. **Open** the sample pattern: `.xy`, `.xye`, `.txt`, `.dat`, `.csv`, `.asc` or `.xrdml`. For XRDML files, the wavelength recorded by the instrument is used.
3. **Analyse.** Völva then runs these steps:
   - fits the background;
   - detects peaks;
   - fits every peak as a pseudo-Voigt Kα₁ + Kα₂ doublet, fitting overlapping peaks together;
   - removes the instrument profile;
   - runs all six methods.
4. **Check the fits.** Use the *Pattern*, *Peak* (◀ ▶) and *Gallery* tabs. Untick doubtful peaks in the table, and type hkl labels if known; the results update immediately.
5. **Read the *Result* tab.** It gives a headline value, what the methods say together, and the full method table.
6. **File ▸ Export results.** This writes a report, CSV tables, seven figures and a session file (see *Outputs*).

### Command line (batch)

```bat
python -m volva calibrate LaB6_2026-09.xy -o Aeris_Cu_quarter.json --note "Aeris, Cu, 1/4° DS" --figure
python -m volva analyse data\*.xy --profile Aeris_Cu_quarter.json -o results
```

---

## What Völva does, and why

**Instrumental correction is Voigt-based.** Each fitted pseudo-Voigt (FWHM, η) is split into Gaussian (G) and Lorentzian (L) widths using the Thompson–Cox–Hastings relations. The instrument profile is modelled as Rietveld programs do: the total FWHM of the standard follows a Caglioti curve, FWHM² = U tan²θ + V tanθ + W, and η is constant (or linear in 2θ if the trend is significant); instrument G and L at any angle come from TCH decomposition of these smooth curves. Standard and sample peaks are fitted as split pseudo-Voigts, so axial-divergence asymmetry does not bias widths. Sample widths are then obtained correctly: G in quadrature, L linearly. Völva 4.0 subtracted total FWHMs in quadrature, which underestimated sizes by 7–35% for 20–200 nm crystallites.

**Breadths are integral breadths.** With K = 1, the Scherrer equation then gives the volume-weighted column length L_vol, the quantity certified for NIST SRM 1979. *Sphere diameter* (K = 4/3) and a custom K are also available.

**Uncertainties are propagated.** Fit covariances and the calibration scatter are carried through by Monte Carlo, using a fixed seed so identical inputs give identical outputs. Linear methods report regression standard errors, and per-reflection methods report the spread between reflections.

**Calibration error is treated as a systematic.** An error in the instrument profile is common to every reflection, so neither of the statistics above captures it. Völva therefore re-runs every method with the instrument profile shifted by ±1σ, both uniformly and as a tilt in tan θ (which can mimic strain). The headline quotes the resulting systematic uncertainty. A strain is only called resolved if it exceeds 2 × (statistical ⊕ calibration) uncertainty.

**Sample components need to be significant.** A Gaussian or Lorentzian sample component is assigned only if the observed width exceeds the instrument's by more than 2σ. Otherwise the whole excess breadth goes to the other component. This stops noise in the quadrature subtraction from inflating breadths near the resolution limit.

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
| Scherrer | total sample breadth, each reflection | — | all broadening is size |
| Single-line Voigt (de Keijser 1982) | Lorentzian part, each reflection | Gaussian part, each reflection | size is Lorentzian, strain is Gaussian |
| Monshi-Scherrer (2012) | intercept of ln β vs ln(1/cos θ) | slope ≠ 1 signals strain | size-only |
| Williamson-Hall (UDM) | intercept | slope | size and strain breadths add linearly |
| Voigt Williamson-Hall (Langford 1992) | Lorentzian intercept | Gaussian slope | size is Lorentzian, strain is Gaussian |
| Halder-Wagner (1966) | slope | intercept | parabolic Voigt approximation |

The *Result* tab turns these contrasts into explicit statements:

- **Strain present?** Völva asks whether W-H, Voigt W-H, H-W and single-line strains are both significant (above 2σ) and material (at least 15% of the high-angle breadth). It also checks whether Scherrer sizes fall with sin θ. If strain is present, Scherrer sizes are lower bounds and the headline switches to a size–strain method.
- **Anisotropic?** Per-reflection sizes are tested for spread beyond their uncertainties (χ² test, CV > 20%). If they differ, a single size is not representative; Völva says so and reports the range. Add hkl labels to see which directions are short.
- **Monshi slope** consistent with 1 means size-only broadening.
- **Profile shape.** A mostly Lorentzian sample profile points to size broadening; a mostly Gaussian one points to strain or defects.
- **Agreement.** The spread of sizes across methods is shown against their differing assumptions.

The size–strain plot (SSP) is not reported separately. In its dimensionally consistent form it is algebraically identical to Halder-Wagner (Pal, *J. Alloys Compd.* **1048**, 185324, 2025); the unit test `test_ssp_dimensionally_consistent_form_equals_hw` checks this.

---

## Outputs (File ▸ Export results)

| File | Contents |
|---|---|
| `<name>_report.txt` | headline, findings, method table, per-reflection table, method notes, instrument profile, full log |
| `<name>_peaks.csv` | per reflection: position, FWHM, η, G/L (observed, instrument, sample), breadths, ratio, Scherrer size or bound, flags |
| `<name>_summary.csv` | per method: size ± error or bound, strain ± error or limit, R², n, warnings |
| `<name>_session.json` | complete record: software and library versions, SHA-256 of sample and standard files, all settings, profile parameters, results, log |
| `<name>_*.png` | pattern, sizes by reflection, comparison, W-H, Voigt W-H, H-W, instrument (300 dpi) |

Every CSV starts with `#` header lines giving the version, run time, file checksum, instrument profile, wavelength and K.

---

## Settings

| Setting | Default | Notes |
|---|---|---|
| Radiation | Cu Kα | *From file* for XRDML; *Custom* for other sources |
| Fit Kα₂ | on | Turn off only for monochromatic (Kα₁) data. It must match the calibration; Völva warns if not |
| Size reported as | L_vol (K = 1) | Sphere diameter (K = 4/3), or custom K |
| Background half-width | automatic | SNIP clipping width in **degrees**, set from peak widths. Völva warns if it is smaller than 1.5 × the widest FWHM |
| Detection threshold | 5 σ | Net peak area above a local baseline, in standard deviations of the local counting noise. Features narrower than the instrument allows are ignored |
| Minimum fit R² | 0.90 | Group fit quality needed to accept a peak. A weak peak with lower R² is still accepted if its residuals are within counting noise (χ²ν ≤ 1.5). Every accepted peak also needs a fitted area of at least 5σ and a width no narrower than 0.7 × the instrument FWHM |
| Resolution limit | 1.2 | β_obs/β_inst below this flags a reflection as near the resolution limit. If most reflections are, only a lower bound is reported |

---

## Tests

`python -m pytest` runs 42 tests (about one minute):

- **Profile mathematics:** FWHM = 2σ for lmfit's pseudo-Voigt; TCH round-trip; Voigt integral breadth against numerical integration.
- **File formats:** headers, comma decimals, delimiters, XRDML with wavelength.
- **Method algebra:** exact recovery from noise-free breadths; SSP ≡ H-W.
- **End-to-end recovery** from synthetic patterns (true Voigt convolution, Kα doublet, Poisson noise):

| Case | Result |
|---|---|
| Strain-free, 10–120 nm | Scherrer L_vol within −0.4% to +2.1% |
| e = 0.05–0.15%, 30–120 nm | size within −9% to 0%; strain within −10% to +2% (Voigt W-H) |
| Strain-free, 340–475 nm (near the resolution limit) | within −2% to +1%, with a calibration systematic of ± 30–70 nm and a measurable-range warning |
| ≥ about 600 nm on this synthetic instrument | reported only as a 2σ lower bound |
| Step size 0.01° / 0.02° / 0.03° (5 nm crystallites) | results agree within 6% |
| 1–5 µm crystallites | reported as "> X nm" (always below the true size), never as a number |
| Other checks | anisotropy detected when real, not flagged for noisy data; missing profile flagged; wrong wavelength caught by the LaB₆ position check; Kα₂ / radiation mismatch blocks the result; single-reflection and odd-setting edge cases; identical input gives identical output |

Two caveats apply. The synthetic data use the same Voigt model that Völva assumes, so real validation (`VALIDATION.md`) is still required. Small strain in very small crystallites (10 nm, e ≤ 0.15%) is often not resolved; this is the detection limit of integral-breadth methods, and Völva reports it as mixed or unresolved rather than guessing.

---

## Limitations

- Integral-breadth methods assume symmetric profiles. They cannot separate size distributions from shape, or treat stacking faults explicitly. For those, use whole-powder-pattern modelling (e.g. WPPM, or TOPAS / GSAS-II with microstructure models).
- Strain is isotropic in all linear methods. Anisotropic strain shows up as scatter and is reported as such, not modelled.
- Peak asymmetry (axial divergence at low angle) is not modelled; low-angle reflections below about 20° 2θ deserve extra scrutiny.
- Binary formats (.raw, .brml) are not read. Export as .xy or .xrdml, or use Bifrost.

## References

- Caglioti, Paoletti & Ricci (1958) *Nucl. Instrum.* 3, 223.
- de Keijser, Langford, Mittemeijer & Vogels (1982) *J. Appl. Cryst.* 15, 308.
- Halder & Wagner (1966) *Acta Cryst.* 20, 312.
- Langford (1992) *NIST Spec. Publ.* 846, 110.
- Monshi, Foroughi & Monshi (2012) *World J. Nano Sci. Eng.* 2, 154.
- Pal (2025) *J. Alloys Compd.* 1048, 185324.
- Ryan et al. (1988) *Nucl. Instrum. Methods B* 34, 396 (SNIP).
- Stokes & Wilson (1942) *Proc. Camb. Phil. Soc.* 38, 313.
- Thompson, Cox & Hastings (1987) *J. Appl. Cryst.* 20, 79.
- Williamson & Hall (1953) *Acta Metall.* 1, 22.
- NIST SRM 660c (LaB₆) and SRM 1979 (ZnO) certificates.

University of Sheffield, N. Reeves-McLaren.
