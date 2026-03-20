# scripts/

Standalone analysis scripts used to produce the results reported in the
manuscript.  These are *not* part of the installable R package — they are
the driver scripts that call package functions on the actual data.

## Files

| File | Description |
|------|-------------|
| `fapa_main.R` | Full analysis driver: loads data, runs all three verification stages, computes BCa CIs, writes output CSVs and plots |

## Usage

1. Install the FAPA package:
   ```r
   remotes::install_github("sekangakim/FAPA")
   ```
2. Place `Calibration.csv` in your working directory (see `data-raw/README.md`).
3. Edit the `USER INPUTS` section at the top of `fapa_main.R` to match your
   paths and preferences.
4. Source the script:
   ```r
   source("scripts/fapa_main.R")
   ```

## Output files produced

| File | Contents |
|------|----------|
| `fapa_CoreProfile1.csv` | BCa CI table for Core Profile 1 |
| `fapa_CoreProfile2.csv` | BCa CI table for Core Profile 2 |
| `fapa_Stage1_PA.csv` | Parallel analysis results |
| `fapa_Stage2_Procrustes.csv` | Principal angle results |
| `fapa_Stage3_TuckerCC.csv` | Tucker CC results |
| `fapa_PersonWeights.csv` | Person-level reconstruction weights |
| `fapa_PersonWeights_Bootstrap.csv` | Bootstrap CIs for selected persons |
