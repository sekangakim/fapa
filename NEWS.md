# FAPA 0.1.0 (2026-03-20)

## Initial release

* `fapa_core()`: SVD-based core profile estimation with sign normalisation
  and variance decomposition.
* `load_and_ipsatize()`: CSV loading and row-centering.
* Three-stage bootstrap verification framework:
  - Stage 1 `fapa_pa()`: variance-matched Horn's parallel analysis.
  - Stage 2 `fapa_procrustes()`: subspace stability via Procrustes
    principal angles.
  - Stage 3 `fapa_tucker()`: profile replicability via Tucker's
    congruence coefficients.
* `fapa_bca()`: BCa bootstrap confidence intervals for core-profile
  coordinates via `boot::boot.ci()`.
* `fapa_person()`: person-level reconstruction, pattern weights, and
  optional bootstrap CIs for selected participants.
* Formatted print methods: `print_pa()`, `print_procrustes()`,
  `print_tucker()`.
* Plotting functions: `plot_fapa_core()`, `plot_pa_scree()`,
  `plot_person_match()`, `plot_principal_angles()`, `plot_tucker_cc()`.
* CSV output helpers: `write_fapa_results()`, `write_verification()`.
