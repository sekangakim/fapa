# CRAN submission comments — fapa 0.1.0

## R CMD check results
0 errors | 0 warnings | 0 notes

## Test environments
- macOS Tahoe 26.3.1, R 4.5.2 (aarch64-apple-darwin20) — local
- Windows (R release) — win-builder
- Ubuntu (R release and devel) — rhub

## Submission notes
- First CRAN submission.
- Synthetic dataset included; no real clinical data distributed with the package.

## Notes for CRAN reviewers
- The `data/fapa_simdata.rda` dataset is a fully synthetic simulation;
  no real clinical data are included in the package.
- The `Calibration.csv` file referenced in the vignette is confidential
  clinical data and is intentionally absent from the package.  The vignette
  uses `eval = FALSE` for all data-dependent chunks and can be built without
  the raw data.
- `boot` is listed as a hard `Imports` dependency because `fapa_bca()` calls
  `boot::boot()` and `boot::boot.ci()` on every run. `boot` is a recommended
  package that ships with every standard R installation, so no additional
  installation burden is imposed on users.

## Downstream dependencies
None at this time (new submission).
