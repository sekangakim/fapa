# FAPA: Factor Analytic Profile Analysis of Ipsatized Data

<!-- badges: start -->
[![R-CMD-check](https://github.com/sekangakim/FAPA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sekangakim/FAPA/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**FAPA** implements a metric inferential framework for pattern detection and
person-level reconstruction in multivariate profile data.  After removing
overall profile elevation via row-centering (ipsatization), FAPA applies
singular value decomposition (SVD) to recover shared *core profiles* and
individual pattern weights, and provides a unified three-stage bootstrap
verification framework.

---

## Motivation

Many existing approaches to profile analysis — including PCA applied to
ipsatized data — lack an integrated, workflow-level treatment of
dimensionality selection, subspace stability, and profile replicability.
FAPA addresses this by combining well-validated tools (Horn's parallel
analysis, Procrustes principal angles, Tucker's congruence coefficients, BCa
bootstrap CIs via the `boot` package) into a coherent, reproducible pipeline
designed specifically around the ipsatized structure of the data.

---

## Three-stage verification framework

| Stage | Method | Criterion |
|-------|--------|-----------|
| 1 | Horn's parallel analysis (variance-matched) | Retain components where σ²_k exceeds the (1−α) quantile of a variance-matched null |
| 2 | Procrustes / principal angles | All K principal angles < 30° across bootstrap resamples |
| 3 | Tucker's congruence coefficients | Mean CC ≥ 0.85 (fair) or ≥ 0.95 (high) per Lorenzo-Seva & ten Berge (2006) |

---

## Installation

```r
# Install from GitHub (development version):
# install.packages("remotes")
remotes::install_github("sekangakim/FAPA")
```

---

## Quick start

```r
library(FAPA)

# 1. Load and ipsatize
dat    <- load_and_ipsatize("Calibration.csv", col_labels = paste0("V", 1:22))
Xtilde <- dat$ipsatized

# 2. Stage 1: parallel analysis
pa  <- fapa_pa(Xtilde, B = 2000, seed = 1)
print_pa(pa)
K   <- pa$n_retain

# 3. Core solution
fit <- fapa_core(Xtilde, K = K)

# 4. Stages 2 & 3: stability and replicability
pr  <- fapa_procrustes(Xtilde, K = K + 3, B = 2000, seed = 1)
tc  <- fapa_tucker(Xtilde,     K = K + 3, B = 2000, seed = 1)
print_procrustes(pr, K_pa = K)
print_tucker(tc, cc_thresh = 0.85, K_pa = K)

# 5. BCa confidence intervals
bca <- fapa_bca(Xtilde, K = K, B = 2000, seed = 1)
plot_fapa_core(bca, i = 1, split_at = 11)

# 6. Person-level reconstruction
prs <- fapa_person(Xtilde, fit, participants = 1:5, seed = 1)
cat(sprintf("Mean R² = %.4f\n", prs$R2_mean))
```

See `vignette("fapa_example", package = "FAPA")` for a fully documented
walkthrough using EDI-2 eating-disorder symptom data.

---

## Key functions

| Function | Description |
|----------|-------------|
| `load_and_ipsatize()` | Load CSV and row-centre |
| `fapa_core()` | SVD-based core profile estimation |
| `fapa_pa()` | Stage 1: variance-matched parallel analysis |
| `fapa_procrustes()` | Stage 2: subspace stability via principal angles |
| `fapa_tucker()` | Stage 3: Tucker's congruence coefficients |
| `fapa_bca()` | BCa bootstrap CIs for core profiles |
| `fapa_person()` | Person-level reconstruction and weights |
| `print_pa/procrustes/tucker()` | Formatted console summaries |
| `plot_fapa_core()` | Core profile with BCa CI bands |
| `plot_pa_scree()` | PA scree plot |
| `plot_principal_angles()` | Stage 2 angle distributions |
| `plot_tucker_cc()` | Stage 3 CC distributions |
| `write_fapa_results()` | Write core-profile CSVs |
| `write_verification()` | Write three-stage verification CSVs |

---

## Citation

If you use FAPA in your research, please cite the accompanying manuscript:

> Kim, S.-K. (in preparation). *Factor Analytic Profile Analysis of
> Ipsatized Data (FAPA): An inferential framework for pattern detection
> and person reconstruction*.

and the software:

> Kim, S.-K. (2026). *FAPA: Factor Analytic Profile Analysis of Ipsatized
> Data*. R package version 0.1.0.
> https://github.com/sekangakim/FAPA

---

## References

- Davison, A. C., & Hinkley, D. V. (1997). *Bootstrap Methods and Their
  Application*. Cambridge University Press.
- Horn, J. L. (1965). A rationale and test for the number of factors in
  factor analysis. *Psychometrika, 30*(2), 179–185.
- Kim, S.-K. (2023). Factorization of person response profiles to identify
  summative profiles carrying central response patterns. *Psychological
  Methods*. https://doi.org/10.1037/met0000568
- Lorenzo-Seva, U., & ten Berge, J. M. F. (2006). Tucker's congruence
  coefficient as a meaningful index of factor similarity. *Methodology,
  2*(2), 57–64.
- Ten Berge, J. M. F. (1999). A legitimate case of component analysis of
  ipsative measures. *Multivariate Behavioral Research, 34*(1), 89–102.

---

## License

MIT © Se-Kang Kim
