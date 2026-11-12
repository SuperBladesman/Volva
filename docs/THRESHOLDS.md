# Thresholds in the code

The acceptance defaults that users can change are in `volva/constants.py`, each marked as sourced or unsourced, and tabulated in the README. The thresholds below are still written into the code where they are used. None of them can be changed without editing the code, and none has a source beyond the statistical conventions noted. Each is marked the same way.

This list was taken from 4.1.2. Moving the ones that decide results into `constants.py` is an open issue.

## Detection and fitting

| Value | Where | What it decides | Source |
|---|---|---|---|
| 0.5 × the minimum width | `peaks.py`, `detect_peaks` | A detected feature narrower than half the minimum fitted width (below) is ignored before fitting | Unsourced (4.1.2) |
| 0.7 × instrument FWHM, never below 1.5 scan steps | `pipeline.py`, `run_analysis` | The narrowest a fitted sample peak may be; 1.5 steps alone for standards and without a profile | 0.7 is `MIN_FWHM_FRACTION_OF_INSTRUMENT`; 1.5 steps unsourced |
| 0.5 × the larger FWHM | `peaks.py`, `_merge_split` | Two maxima closer than this are one peak with a noisy top | Unsourced (4.1.2) |
| 0.35 × the Kα₁–Kα₂ separation (at least 2 steps), and 0.85 × the height | `peaks.py`, `_drop_ka2` | A weaker detection at the Kα₂ position of a peak is that peak's Kα₂ line, not a reflection | Unsourced (4.1.0, audit item M1) |
| η ≤ 0.005 or ≥ 0.995 | `peaks.py`, `_fit_once`; also `instrument.py` and `diagnostics.py` | η is at a bound, and so not determined | Unsourced |
| Width within 1 % of a fitting bound | `peaks.py`, `_fit_once` | The fit is rejected (`width_at_bound`) | Unsourced |
| Another peak within 1.5 × FWHM | `peaks.py`, `_fit_once` | The peak is flagged as overlapped | Unsourced |
| 0.05 | `peaks.py`, `_fit_once`; `correction.py`, `correct_peaks` | The uncertainty given to η when it is fixed at a bound | Unsourced |
| 5 % of the FWHM and 0.08 in η | `correction.py`, `correct_peaks` | The fit uncertainties assumed when a fit returns no covariance (flagged `assumed_fit_uncertainty`) | Unsourced |
| 0.25 × the smaller FWHM | `peaks.py`, `fit_all` | Two fits closer than this are one peak detected twice | Unsourced |
| 1.4826 × the median absolute deviation | `background.py`, `local_noise`; `instrument.py`, `calibrate_from_standard` | Converts a median absolute deviation into a standard deviation | Sourced: the consistency constant for normally distributed data (Rousseeuw & Croux, *J. Am. Stat. Assoc.* 88, 1273, 1993) |

## Instrument calibration

| Value | Where | What it decides | Source |
|---|---|---|---|
| 4 × the scaled median absolute deviation, up to three rounds | `instrument.py`, `calibrate_from_standard` | A standard peak this far from the Caglioti curve is left out of the profile | Unsourced |
| Slope > 2 standard errors | `instrument.py`, `calibrate_from_standard` | A linear trend in η is kept; otherwise η is constant | Unsourced (statistical convention, about 95 %) |
| 0.02 | `instrument.py`, `calibrate_from_standard` | The smallest uncertainty allowed for the η of a standard peak | Unsourced |
| FWHM scatter > 10 % of the median width; median FWHM < 3 steps; η scatter > 0.08, or two or more peaks at a bound | `instrument.py`, `calibrate_from_standard` | Calibration warnings. The advice given, a step of at most FWHM/5, is also unsourced | Unsourced |
| Peaks within 0.3° of a reflection; at least 60 % matched; mean offset < 0.1° | `instrument.py`, `_position_check` | Whether the standard, the wavelength and the zero error are consistent | Unsourced |

## Instrumental correction

| Value | Where | What it decides | Source |
|---|---|---|---|
| Excess over the instrument > 2σ (fit and calibration combined) | `correction.py`, `correct_peaks` | A Gaussian or Lorentzian sample component is assigned; otherwise the whole excess breadth goes to the other component | Unsourced (statistical convention) |
| Sample breadth < 2σ | `correction.py`, `correct_peaks` | The reflection is resolution-limited: it gives only a lower bound on size (from β + 2σ) and is left out of means | Unsourced (statistical convention) |

## Methods

| Value | Where | What it decides | Source |
|---|---|---|---|
| Component breadth > 2σ | `methods.py`, `single_line` | A reflection's Lorentzian part gives a size, or its Gaussian part a strain | Unsourced (statistical convention) |
| Most reflections strained, mean > 2 standard errors, and CV < 50 % | `methods.py`, `single_line` | Single-line strain is significant | Unsourced |
| \|t\| > 2 for slope − 1 | `methods.py`, `monshi` | The Monshi slope differs from 1, so broadening is not size-only and the Monshi size is not used | Unsourced (statistical convention) |
| \|slope\| ≤ 2σ | `methods.py`, `_strain_from_slope` | A linear-method strain is not significant and is reported as an upper limit (slope + 2σ) | Unsourced (statistical convention) |

## Interpretation

| Value | Where | What it decides | Source |
|---|---|---|---|
| 4e tan θ ≥ 15 % of the sample breadth, over the three highest-angle reflections | `diagnostics.py`, `_interpret` | A significant strain is also material; otherwise it is treated as negligible | Unsourced |
| Strain > 2 × (statistical ⊕ calibration uncertainty) | `pipeline.py`, `_systematic` | A strain survives the calibration systematic; otherwise it is an upper limit | Unsourced (statistical convention) |
| r < −0.5 and p < 0.05 | `diagnostics.py`, `_interpret` | Scherrer sizes fall with angle, which counts as evidence of strain | Unsourced |
| CV > 20 % and χ² p < 0.01 over at least three reflections; with four or more, it must survive removal of the most deviant one | `diagnostics.py`, `_interpret` | The sample is reported as anisotropic | Unsourced (4.1.0; the single-reflection check 4.1.1) |
| Fewer than two reflections, or fewer than half of them, with η free and σ(η) < 0.25 | `diagnostics.py`, `_interpret` | The Gaussian/Lorentzian split is not determined, so single-line Voigt and Voigt Williamson–Hall are not used | Unsourced (4.1.2) |
| Mean Lorentzian fraction > 0.75, or < 0.35 | `diagnostics.py`, `_interpret` | The profile is described as mostly Lorentzian (size) or mostly Gaussian (strain) | Unsourced |
| Scatter > 15 % about the line, with at least four points | `diagnostics.py`, `_interpret` | Williamson–Hall or Halder–Wagner is not linear enough for the isotropic model | Unsourced |
| Spread < 15 % of the median, with at least three methods | `diagnostics.py`, `_interpret` | The methods are said to agree | Unsourced |
| Scherrer size > 0.5 × the largest resolvable size | `diagnostics.py`, `_interpret` | Measurable-range warning | Unsourced |
| sin θ spans < 0.25 | `diagnostics.py`, `_interpret` | Angular-range warning for the linear methods | Unsourced |
| Relative difference in λ(Kα₁) > 2 × 10⁻⁴ | `pipeline.py`, `_compatibility` | The profile and the sample were measured at different wavelengths; no values are reported | Unsourced |
| SNIP half-width < 1.5 × the widest FWHM | `pipeline.py`, `run_analysis` | Warning that the background may cut into peak tails | Unsourced |
