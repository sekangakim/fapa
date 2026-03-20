# data-raw/

Place raw data preparation scripts and original data files here.
**This directory is excluded from the built R package** (see `.Rbuildignore`).

## Data files used in the manuscript

| File | Description | Status |
|------|-------------|--------|
| `Calibration.csv` | EDI-2 pre/post data, n = 2,599 persons × 22 variables | Confidential — not distributed |
| `Validation.csv`  | Hold-out sample for split-sample replication | Confidential — not distributed |

## Why data are not included

The EDI-2 data were provided by Bunmi Olatunji (Vanderbilt University)
for research purposes and cannot be redistributed.  To reproduce the
manuscript analyses, researchers should use their own ipsatized
person-by-variable data and adapt `scripts/fapa_main.R` accordingly.

## Simulated example data

A synthetic dataset approximating the structure of the calibration sample
is provided as the package dataset `fapa_simdata`.  To regenerate it:

```r
source("data-raw/simulate_fapa_data.R")
usethis::use_data(fapa_simdata, overwrite = TRUE)
```
