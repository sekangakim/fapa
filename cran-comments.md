# CRAN submission comments — FAPA 0.1.1

## Resubmission
This is a resubmission in response to reviewer feedback (Konstanze Lauseker).
The following changes were made:

- Added single quotes around software and package names ('FAPA', 'SVD', 'boot')
  in the Description field of DESCRIPTION, as requested.
- Added a doi reference for Davison and Hinkley (1997) in DESCRIPTION:
  <doi:10.1017/CBO9780511802843>.
- Replaced \dontrun{} with \donttest{} in the load_and_ipsatize() example,
  as the example requires an external file (Calibration.csv) that is
  intentionally absent from the package due to data confidentiality.
- Added 'utils' to the Imports field in DESCRIPTION to formally declare the
  dependency on utils::read.csv() and utils::write.csv().

## R CMD check results
0 errors | 0 warnings | 0 notes

## Test environments
- macOS Tahoe 26.3.1, R 4.5.2 (aarch64-apple-darwin20) — local
- Windows Server 2022, R-devel 4.6.0 alpha (win-builder): 0 errors,
  0 warnings, 1 note (new submission only)

## Notes for CRAN reviewers
- The `data/fapa_simdata.rda` dataset is a fully synthetic simulation;
  no real clinical data are included in the package.
- The `Calibration.csv` file referenced in the vignette is confidential
  clinical data and is intentionally absent from the package. The vignette
  uses `eval = FALSE` for all data-dependent chunks and can be built without
  the raw data. The load_and_ipsatize() example is wrapped in \donttest{}
  for this reason.
- 'boot' is listed as a hard Imports dependency because fapa_bca() calls
  boot::boot() and boot::boot.ci() on every run. 'boot' is a recommended
  package that ships with every standard R installation, so no additional
  installation burden is imposed on users.

## Downstream dependencies
None at this time (new submission).
