## data-raw/simulate_fapa_data.R
## Generates fapa_simdata: a synthetic EDI-2-like dataset for package examples.
## Approximates the calibration sample (n = 2599, p = 22) in structure but
## contains no real clinical records.
##
## Run once to regenerate, then call usethis::use_data(fapa_simdata).

set.seed(20260320)

n_persons  <- 500   # smaller than real data but sufficient for examples
n_vars     <- 22    # 11 pre + 11 post EDI-2 subscales

## EDI-2 subscale tags
edi_tags      <- c("Dt","Bu","Bd","In","Pf","Id","Ia","Mf","As","Ir","Si")
before_labels <- paste0("Before_",  1:11, "_", edi_tags)
after_labels  <- paste0("After_",  12:22, "_", edi_tags)

## Approximate marginal means and SDs from the calibration sample
## (Before subscales generally higher than After for a treatment sample)
mu_before <- c(18, 12, 22, 14, 16, 10, 12, 8,  9,  11, 10)
mu_after  <- c(14,  9, 17, 11, 14,  8, 10, 7,  7,   9,  8)
sd_vals   <- rep(6, n_vars)

## Latent structure: 2 components (mirrors K_pa = 2 finding)
## Component 1: normative gradient (high Dt, Bd; low Mf, As)
## Component 2: post-treatment change contrast
V1 <- c( 1.2,  0.8,  1.5,  0.6,  0.9,  0.4,  0.5, -0.3, -0.2,  0.3,  0.2,
          0.9,  0.6,  1.2,  0.5,  0.7,  0.3,  0.4, -0.2, -0.1,  0.2,  0.1)
V2 <- c( 0.3,  0.2,  0.4,  0.1,  0.2,  0.0,  0.1,  0.1,  0.0,  0.1,  0.1,
         -0.8, -0.5, -1.0, -0.4, -0.6, -0.2, -0.3,  0.1,  0.1, -0.2, -0.1)

## Generate data with latent structure plus noise
scores_latent <- outer(rnorm(n_persons), V1) +
                 outer(rnorm(n_persons), V2)
noise         <- matrix(rnorm(n_persons * n_vars, sd = 2), n_persons, n_vars)
raw_scores    <- scores_latent + noise

## Shift to approximate marginal means and clip to plausible item range
mu_vec     <- c(mu_before, mu_after)
raw_scores <- sweep(raw_scores, 2, colMeans(raw_scores), "-")
raw_scores <- sweep(raw_scores, 2, mu_vec, "+")
raw_scores <- round(pmax(pmin(raw_scores, 40), 0))   # EDI-2 range 0-40

fapa_simdata <- as.data.frame(raw_scores)
colnames(fapa_simdata) <- c(before_labels, after_labels)

## Save as package data
usethis::use_data(fapa_simdata, overwrite = TRUE)
message(sprintf("fapa_simdata saved: %d persons x %d variables",
                nrow(fapa_simdata), ncol(fapa_simdata)))
