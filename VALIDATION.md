# Völva 4.1: validation protocol

**Purpose:** show, on this lab's instruments, that Völva reports crystallite size within agreed limits and never reports a size for crystallites beyond instrumental resolution. Complete this protocol before the tool is used for reported results.

## 1. Before measuring: settle the convention

Völva reports **L_vol**, the volume-weighted mean column length (integral breadth, K = 1). NIST SRM 1979 certifies L_vol for each reflection, so no conversion is needed. Do not compare with FWHM-based sizes or with K = 0.9/0.94 values.

## 2. Materials

| Material | Role | Check |
|---|---|---|
| NIST SRM 660c LaB₆ | Instrument profile | FWHM scatter about the Caglioti curve ≤ 5% of the typical width; no standard peak more than 10% from the model; no η-noise warning; position check passes. Scan with a step ≤ FWHM/5 (≈ 0.015° on a typical Bragg–Brentano instrument) |
| NIST SRM 1979 ZnO, 15 nm material | Small-size accuracy | Per-line L_vol against certificate (certified 14.5–28.6 nm across 22 lines) |
| NIST SRM 1979 ZnO, 60 nm material | Larger-size accuracy | Per-line L_vol against certificate (certified 66.7–112.4 nm) |
| CeO₂ round-robin sample (Balzar et al., *J. Appl. Cryst.* 37, 911, 2004), if available | Near-isotropic check | Size and strain against the published consensus |
| Well-crystallised LaB₆ or Si run as a *sample* | Negative control | Must be reported as "> X nm (resolution-limited)" |
| Two or three routine lab samples | Real-world agreement | Compare with whole-pattern fitting (GSAS-II / HighScore) |

The SRM 1979 certification used SRM 660b for its instrument profile. Record which LaB₆ lot you use.

## 3. Measurements

Use the lab's routine optics and scan programme; the instrument profile must match what users will run. For each material:

1. Take three fresh mounts, by two operators.
2. Measure at step sizes of 0.01°, 0.02° and 0.03° (same counting statistics per degree where possible).
3. Analyse with default settings, then again after unticking any flagged peaks. Record both.

Save each run's export bundle; the `_session.json` file is the record.

## 4. Suggested acceptance criteria (to be agreed with the lab manager)

| Criterion | Limit |
|---|---|
| 60 nm ZnO, per-line L_vol | within ±10% of the certificate |
| 15 nm ZnO, per-line L_vol | within ±15% of the certificate |
| Repeatability (mounts, operators) | CV < 5% |
| Step-size dependence | < 5% change in L_vol |
| Negative control | flagged resolution-limited in 100% of runs |
| Upper reporting limit | largest L_vol the lab will quote; set from the 60 nm results and the resolution threshold |

## 5. Results log (fill in)

| Date | Material | Mount / operator | Step (°) | Profile file | Völva version | Headline | Per-line within limit (n/N) | Pass? |
|---|---|---|---|---|---|---|---|---|
| | | | | | 4.1.2 | | | |

## 6. Re-validation triggers

Re-run the SRM 1979 check after any of the following:

- a change of optics, detector or tube;
- a new instrument profile made with a different standard lot;
- any Völva release that changes `profiles.py`, `peaks.py`, `correction.py` or `methods.py`.
