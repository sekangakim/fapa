# CRAN submission comments — FAPA 0.1.0

## R CMD check results

0 errors | 0 warnings | 0 notes

Tested on:
- macOS (R release)
- Windows (R release)
- Ubuntu (R release and devel)

## Notes for CRAN reviewers

- The `data/fapa_simdata.rda` dataset is a fully synthetic simulation;
  no real clinical data are included in the package.
- The `Calibration.csv` file referenced in the vignette is confidential
  clinical data and is intentionally absent from the package.  The vignette
  uses `eval = FALSE` for all data-dependent chunks and can be built without
  the raw data.
- `boot` is listed as a hard `Imports` dependency because `fapa_bca()` calls
  `boot::boot()` and `boot::boot.ci()` on every run.  `boot` ships with base R
  so no new installation burden is imposed on users.

## Downstream dependencies

None at this time (new submission).
