## =============================================================================
## fapa_main.R
## Factor Analytic Profile Analysis of Ipsatized Data (FAPA)
## Se-Kang Kim
##
## Three-stage bootstrap verification framework:
##   Stage 1 – Horn's parallel analysis  (dimensionality)
##   Stage 2 – Procrustes / principal angles (subspace stability)
##   Stage 3 – Tucker's congruence coefficient (profile replicability)
##
## All verification stages use B = 2000 bootstrap replicates.
## =============================================================================

suppressPackageStartupMessages({
  require(stats)
  require(graphics)
  require(boot)          # boot() + boot.ci(type="bca")  — Davison & Hinkley (1997)
})

## =============================================================================
## USER INPUTS  (edit this section only)
## =============================================================================

setwd("/Users/sekangkim/Documents/J Sub/FAPA/")

data_path    <- "Calibration.csv"   # person × variable data (no header row assumed)
out_prefix   <- "fapa"              # base name for CSV outputs
seed         <- 1
B_boot       <- 2000                # bootstrap replicates (all three stages)
alpha        <- 0.05                # two-tailed alpha for CIs and PA
angle_thresh <- 30                  # principal-angle stability threshold (degrees): all angles must be < this
cc_thresh    <- 0.85                # Tucker CC lower bound for "acceptable" replication
K_extra      <- 3                   # extra dimensions beyond K_pa to include in Stages 2 & 3 for contrast
participants <- 1:5                 # IDs for individual-level bootstrap inference
direction    <- NULL                # sign corrections (NULL = auto by first positive mean)

## EDI-2 subscale labels (22 variables: 11 Before + 11 After)
edi_tags <- c("Dt","Bu","Bd","In","Pf","Id","Ia","Mf","As","Ir","Si")
before_labels <- paste0("Before_",  1:11, "_", edi_tags)
after_labels  <- paste0("After_",  12:22, "_", edi_tags)


## =============================================================================
## 1.  DATA LOADING AND IPSATIZATION
## =============================================================================

load_and_ipsatize <- function(path, col_labels) {
  dat <- read.csv(path, header = TRUE, check.names = FALSE)
  if (ncol(dat) != length(col_labels))
    stop(sprintf("Expected %d columns, got %d.", length(col_labels), ncol(dat)))
  colnames(dat) <- col_labels
  row_means <- rowMeans(dat, na.rm = TRUE)
  Xtilde    <- as.matrix(sweep(dat, 1, row_means, FUN = "-"))
  list(raw = dat, ipsatized = Xtilde, row_means = row_means,
       varnames = col_labels)
}


## =============================================================================
## 2.  CORE FAPA ESTIMATION
##     SVD of ipsatized matrix → core profiles, person weights, variance
## =============================================================================

fapa_core <- function(Xtilde, K, direction = NULL) {
  ## Thin SVD
  sv  <- svd(Xtilde, nu = K, nv = K)
  U   <- sv$u                           # P × K  person-weight matrix
  S   <- sv$d[seq_len(K)]               # K singular values
  V   <- sv$v                           # I × K  right singular vectors

  ## Core-profile (scale) matrix  X = V Σ
  X   <- V %*% diag(S, nrow = K)        # I × K

  ## Optional sign normalisation: flip so largest-abs coordinate is positive
  if (is.null(direction)) {
    direction <- sign(apply(X, 2, function(x) x[which.max(abs(x))]))
    direction[direction == 0] <- 1
  }
  X <- X %*% diag(direction, nrow = K)
  U <- U %*% diag(direction, nrow = K)
  V <- V %*% diag(direction, nrow = K)

  ## Variance accounting
  total_var <- sum(Xtilde^2)
  var_k     <- S^2
  prop_var  <- var_k / total_var
  cum_var   <- cumsum(prop_var)

  ## Person–core correlations (normalized)
  prof_norms <- sqrt(rowSums(Xtilde^2))
  person_cor <- sweep(U * rep(S, each = nrow(U)), 1,
                       pmax(prof_norms, .Machine$double.eps), "/")

  list(U = U, S = S, V = V, X = X,
       total_var = total_var, var_k = var_k,
       prop_var  = prop_var,  cum_var = cum_var,
       person_cor = person_cor,
       direction  = direction, K = K)
}


## =============================================================================
## STAGE 1 – HORN'S PARALLEL ANALYSIS (variance-matched bootstrap version)
##
## For each of B random datasets of identical dimensions and row-centering,
## the random matrix is rescaled to the same Frobenius norm as Xtilde before
## SVD.  This places observed and random σ² on the same variance scale,
## preventing the trivial outcome where raw-score data always dominates N(0,1)
## random matrices.  Retain components whose observed σ²_k exceeds the
## (1−α) quantile of the matched random distribution.
## =============================================================================

fapa_pa <- function(Xtilde, B = 2000, alpha = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  P     <- nrow(Xtilde);  I <- ncol(Xtilde)
  K_max <- min(P - 1L, I - 1L)

  obs_sv2    <- svd(Xtilde, nu = 0, nv = 0)$d[seq_len(K_max)]^2
  total_var  <- sum(obs_sv2)                 # Frobenius norm² of Xtilde

  rand_sv2 <- matrix(0, B, K_max)
  for (b in seq_len(B)) {
    Xr  <- matrix(rnorm(P * I), P, I)
    Xr  <- Xr - rowMeans(Xr)              # row-centre (ipsatize) random data
    ## Rescale to same total variance as observed ipsatized matrix
    Xr  <- Xr * sqrt(total_var / sum(Xr^2))
    svr <- svd(Xr, nu = 0, nv = 0)$d[seq_len(K_max)]
    rand_sv2[b, ] <- svr^2
  }

  thresh   <- apply(rand_sv2, 2, quantile, probs = 1 - alpha)
  n_retain <- sum(obs_sv2 > thresh)
  if (n_retain == 0L) {
    warning("PA retained 0 components. Check data scaling or reduce alpha.")
    n_retain <- 1L
  }

  prop_obs  <- obs_sv2  / total_var
  prop_rand <- colMeans(rand_sv2) / total_var   # mean random proportion per component

  list(n_retain   = n_retain,
       obs_sv2    = obs_sv2,
       thresh     = thresh,
       prop_obs   = prop_obs,
       prop_rand  = prop_rand,
       rand_sv2   = rand_sv2,
       total_var  = total_var,
       alpha      = alpha,
       B          = B)
}

## Summary printer for Stage 1
print_pa <- function(pa) {
  cat("\n=== STAGE 1: Horn's Parallel Analysis (variance-matched) ===\n")
  cat(sprintf("  Bootstrap replicates : %d\n", pa$B))
  cat(sprintf("  Alpha                : %.3f\n", pa$alpha))
  cat(sprintf("  Components retained  : %d\n\n", pa$n_retain))
  df <- data.frame(
    Component    = seq_along(pa$obs_sv2),
    Obs_sv2      = round(pa$obs_sv2,   3),
    Obs_PropVar  = round(pa$prop_obs,  4),
    Rand_95pct   = round(pa$thresh,    3),
    Rand_PropVar = round(pa$prop_rand, 4),
    Retain       = ifelse(pa$obs_sv2 > pa$thresh, "YES", "no")
  )
  print(df, row.names = FALSE)
}


## =============================================================================
## STAGE 2 – PROCRUSTES / PRINCIPAL ANGLES  (subspace stability)
##
## For each bootstrap replicate, compute the K* right-singular-vector subspace
## of the resampled ipsatized data, then measure the principal angles (in
## degrees) between that subspace and the original V_K*.
##
## Stability criterion: ALL K* principal angles must be < angle_thresh (default 30°).
## A replicate is "stable" only when every angle satisfies this condition,
## confirming that the bootstrap subspace is nearly parallel to the original.
## =============================================================================

## Helper: principal angles between two column spaces (degrees)
principal_angles_deg <- function(A, B) {
  ## A, B: matrices with same number of columns; columns need not be orthonormal
  QA <- qr.Q(qr(A))
  QB <- qr.Q(qr(B))
  cos_ang <- pmin(svd(crossprod(QA, QB), nu = 0, nv = 0)$d, 1)
  acos(cos_ang) * 180 / pi
}

fapa_procrustes <- function(Xtilde, K, B = 2000,
                             angle_thresh = 30, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  P <- nrow(Xtilde);  I <- ncol(Xtilde)
  n <- P

  ## Original subspace (right singular vectors)
  V0 <- svd(Xtilde, nu = 0, nv = K)$v   # I × K

  angles_mat <- matrix(NA, B, K)         # B × K principal angles
  n_stable   <- 0L                       # replicates where ALL angles < threshold

  for (b in seq_len(B)) {
    idx  <- sample(n, n, replace = TRUE)
    Xb   <- Xtilde[idx, , drop = FALSE]
    Vb   <- svd(Xb, nu = 0, nv = K)$v
    ang  <- principal_angles_deg(V0, Vb)
    angles_mat[b, ] <- ang
    if (all(ang < angle_thresh)) n_stable <- n_stable + 1L
  }

  prop_stable <- n_stable / B

  list(angles_mat  = angles_mat,
       angle_mean  = colMeans(angles_mat),
       angle_sd    = apply(angles_mat, 2, sd),
       angle_q025  = apply(angles_mat, 2, quantile, 0.025),
       angle_q975  = apply(angles_mat, 2, quantile, 0.975),
       n_stable    = n_stable,
       prop_stable = prop_stable,
       angle_thresh = angle_thresh,
       K = K, B = B)
}

## Summary printer for Stage 2
print_procrustes <- function(pr, K_pa = NULL) {
  cat("\n=== STAGE 2: Procrustes / Principal Angles ===\n")
  cat(sprintf("  Bootstrap replicates           : %d\n", pr$B))
  cat(sprintf("  Stability threshold            : < %.0f degrees\n", pr$angle_thresh))
  cat(sprintf("  Replicates with ALL angles < threshold: %d (%.1f%%)\n\n",
              pr$n_stable, 100 * pr$prop_stable))
  df <- data.frame(
    Dimension  = seq_len(pr$K),
    Mean_deg   = round(pr$angle_mean, 2),
    SD_deg     = round(pr$angle_sd,   2),
    CI_lo_deg  = round(pr$angle_q025, 2),
    CI_hi_deg  = round(pr$angle_q975, 2),
    Angle_OK   = ifelse(pr$angle_mean < pr$angle_thresh, "YES", "no")
  )
  if (!is.null(K_pa))
    df$PA_Retained <- ifelse(seq_len(pr$K) <= K_pa, "YES", "no")
  print(df, row.names = FALSE)
}


## =============================================================================
## STAGE 3 – TUCKER'S CONGRUENCE COEFFICIENT  (profile replicability)
##
## For each bootstrap replicate, align the bootstrap core profiles X_b to the
## original X via column-wise sign matching (to handle sign ambiguity) and
## compute Tucker's CC per dimension.
##
## Conventional thresholds (Lorenzo-Seva & ten Berge, 2006):
##   CC ≥ 0.95 : high similarity ("factor replication")
##   CC ≥ 0.85 : fair similarity
##   CC < 0.85 : poor similarity
## =============================================================================

tucker_cc <- function(x, y) {
  ## Cosine of the angle between two vectors
  denom <- sqrt(sum(x^2) * sum(y^2))
  if (denom < .Machine$double.eps) return(NA_real_)
  sum(x * y) / denom
}

fapa_tucker <- function(Xtilde, K, B = 2000,
                         cc_thresh = 0.85, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  P <- nrow(Xtilde);  I <- ncol(Xtilde)
  n <- P

  fit0 <- fapa_core(Xtilde, K)
  X0   <- fit0$X                     # I × K  original core profiles

  cc_mat <- matrix(NA, B, K)         # B × K  Tucker CCs

  for (b in seq_len(B)) {
    idx  <- sample(n, n, replace = TRUE)
    Xb   <- Xtilde[idx, , drop = FALSE]
    fitb <- fapa_core(Xb, K)
    Xb_cp <- fitb$X

    for (k in seq_len(K)) {
      ## Resolve sign ambiguity: pick sign that maximises |CC|
      cc_pos <- tucker_cc(X0[, k],  Xb_cp[, k])
      cc_neg <- tucker_cc(X0[, k], -Xb_cp[, k])
      cc_mat[b, k] <- ifelse(abs(cc_pos) >= abs(cc_neg), cc_pos, cc_neg)
    }
  }

  list(cc_mat    = cc_mat,
       cc_mean   = colMeans(cc_mat, na.rm = TRUE),
       cc_sd     = apply(cc_mat, 2, sd, na.rm = TRUE),
       cc_q025   = apply(cc_mat, 2, quantile, 0.025, na.rm = TRUE),
       cc_q975   = apply(cc_mat, 2, quantile, 0.975, na.rm = TRUE),
       cc_thresh = cc_thresh,
       K = K, B = B)
}

## Summary printer for Stage 3
print_tucker <- function(tc, cc_thresh, K_pa = NULL) {
  cat("\n=== STAGE 3: Tucker's Congruence Coefficients ===\n")
  cat(sprintf("  Bootstrap replicates : %d\n", tc$B))
  cat(sprintf("  Acceptability cutoff : CC ≥ %.2f\n\n", cc_thresh))
  df <- data.frame(
    Core_Profile = seq_len(tc$K),
    Mean_CC      = round(tc$cc_mean,  4),
    SD_CC        = round(tc$cc_sd,    4),
    CI_lo        = round(tc$cc_q025,  4),
    CI_hi        = round(tc$cc_q975,  4),
    CC_OK        = ifelse(tc$cc_mean >= cc_thresh, "YES", "no")
  )
  if (!is.null(K_pa))
    df$PA_Retained <- ifelse(seq_len(tc$K) <= K_pa, "YES", "no")
  print(df, row.names = FALSE)
  cat("\n  Reference: Lorenzo-Seva & ten Berge (2006)\n")
  cat("    CC ≥ 0.95 = high; CC ≥ 0.85 = fair; CC < 0.85 = poor\n")
}


## =============================================================================
## 4.  BCa BOOTSTRAP CONFIDENCE INTERVALS FOR CORE PROFILES
##     Uses boot::boot() + boot::boot.ci(type = "bca") — the canonical R
##     implementation (Davison & Hinkley, 1997).  Sign alignment is handled
##     inside the statistic function passed to boot(), using the inner-product
##     rule: flip column k of Xb if <X0[,k], Xb[,k]> < 0.
## =============================================================================

if (!requireNamespace("boot", quietly = TRUE))
  stop("Package 'boot' is required. Install with: install.packages('boot')")

## Internal helper: align columns of Xb to X0 via inner-product sign check.
## Flips column k of Xb if it points away from X0[,k].  This prevents the
## bootstrap distribution of each core-profile coordinate from becoming
## bimodal due to SVD sign ambiguity across resamples.
.align_signs <- function(X0, Xb) {
  for (k in seq_len(ncol(X0)))
    if (sum(X0[, k] * Xb[, k]) < 0) Xb[, k] <- -Xb[, k]
  Xb
}

fapa_bca <- function(Xtilde, K, B = 2000, alpha = 0.05,
                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  P <- nrow(Xtilde);  I <- ncol(Xtilde)

  ## Original solution — defines reference signs for alignment
  fit0 <- fapa_core(Xtilde, K)
  X0   <- fit0$X                       # I × K  original core profiles

  ## Statistic function for boot():
  ##   data    = Xtilde (passed as the 'data' argument)
  ##   indices = bootstrap row indices supplied by boot()
  ##   Returns a vector of length I*K: columns of aligned X stacked
  fapa_stat <- function(data, indices) {
    Xb   <- data[indices, , drop = FALSE]
    fitb <- fapa_core(Xb, K)
    as.vector(.align_signs(X0, fitb$X))   # vectorise column-by-column
  }

  ## Run boot() — handles both bootstrap replicates and jackknife for BCa
  boot_out <- boot::boot(data      = Xtilde,
                         statistic = fapa_stat,
                         R         = B)

  ## Extract BCa intervals for each of the I*K coordinates
  la <- alpha / 2;  ua <- 1 - alpha / 2
  n_stats <- I * K

  bca_lo  <- numeric(n_stats)
  bca_hi  <- numeric(n_stats)
  boot_mn <- numeric(n_stats)
  boot_se <- numeric(n_stats)
  pct_lo  <- numeric(n_stats)
  pct_hi  <- numeric(n_stats)

  for (idx in seq_len(n_stats)) {
    ci <- tryCatch(
      boot::boot.ci(boot_out, conf = 1 - alpha, type = "bca", index = idx),
      error = function(e) NULL
    )
    if (!is.null(ci)) {
      bca_lo[idx] <- ci$bca[4]
      bca_hi[idx] <- ci$bca[5]
    } else {
      ## Fallback to percentile if BCa fails (e.g. constant statistic)
      bca_lo[idx] <- quantile(boot_out$t[, idx], la, na.rm = TRUE)
      bca_hi[idx] <- quantile(boot_out$t[, idx], ua, na.rm = TRUE)
    }
    boot_mn[idx] <- mean(boot_out$t[, idx], na.rm = TRUE)
    boot_se[idx] <- sd(boot_out$t[, idx],   na.rm = TRUE)
    pct_lo[idx]  <- quantile(boot_out$t[, idx], la, na.rm = TRUE)
    pct_hi[idx]  <- quantile(boot_out$t[, idx], ua, na.rm = TRUE)
  }

  ## Reshape results into a list of K data frames (one per core profile)
  ci_list <- vector("list", K)
  for (k in seq_len(K)) {
    idx_k <- seq_len(I) + (k - 1L) * I   # column-k indices in the flat vector
    ci_list[[k]] <- data.frame(
      Ori      = X0[, k],
      Mean     = boot_mn[idx_k],
      SE       = boot_se[idx_k],
      Lower    = pct_lo[idx_k],
      Upper    = pct_hi[idx_k],
      BCaLower = bca_lo[idx_k],
      BCaUpper = bca_hi[idx_k],
      row.names = colnames(Xtilde)
    )
  }

  ## Also expose the full boot object and the raw B × (I*K) matrix
  ## for downstream use (e.g. plots, diagnostics)
  list(ci       = ci_list,
       X0       = X0,
       boot_out = boot_out,
       boot_X   = array(boot_out$t, dim = c(B, I, K)),
       K        = K,
       B        = B,
       alpha    = alpha,
       varnames = colnames(Xtilde))
}


## =============================================================================
## 5.  PERSON-LEVEL RECONSTRUCTION AND WEIGHTS
## =============================================================================

fapa_person <- function(Xtilde, fit, participants = NULL,
                         B_boot = 2000, alpha = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  P <- nrow(Xtilde);  I <- ncol(Xtilde)
  K <- fit$K
  la <- alpha / 2;  ua <- 1 - alpha / 2

  ## Point estimates: reconstruct each person
  U  <- fit$U
  S  <- fit$S
  X  <- fit$X
  row_means <- rowMeans(Xtilde)

  R2_k <- numeric(K)                    # partial R² per dimension
  for (k in seq_len(K)) {
    Xk   <- X[, k, drop = FALSE]
    y_all <- as.matrix(Xtilde)
    R2_k[k] <- 1 - sum((y_all - U[, k, drop = FALSE] %*% t(Xk))^2) /
                    sum(y_all^2)
  }

  ## Full reconstruction R² per person
  Xhat   <- U %*% t(X)
  R2_person <- 1 - rowSums((Xtilde - Xhat)^2) /
                   rowSums(Xtilde^2)

  result <- data.frame(
    Person  = seq_len(P),
    Level   = row_means,
    R2      = round(R2_person, 4),
    matrix(round(U,  4), ncol = K,
           dimnames = list(NULL, paste0("w", seq_len(K)))),
    matrix(round(fit$person_cor, 4), ncol = K,
           dimnames = list(NULL, paste0("rDim", seq_len(K))))
  )

  ## Bootstrap CIs for selected participants only
  result_B <- NULL
  if (!is.null(participants)) {
    n <- P
    boot_w <- array(NA, dim = c(B_boot, length(participants), K))

    for (b in seq_len(B_boot)) {
      idx  <- sample(n, n, replace = TRUE)
      Xb   <- Xtilde[idx, , drop = FALSE]
      fitb <- fapa_core(Xb, K, direction = fit$direction)
      ## Project original participants onto bootstrap core profiles
      for (pi in seq_along(participants)) {
        p  <- participants[pi]
        yp <- as.numeric(Xtilde[p, ])
        xb <- fitb$X
        boot_w[b, pi, ] <- as.numeric(
          solve(crossprod(xb), crossprod(xb, yp)))
      }
    }

    rows <- lapply(seq_along(participants), function(pi) {
      p    <- participants[pi]
      wmat <- boot_w[, pi, ]
      if (K == 1) wmat <- matrix(wmat, ncol = 1)
      c(Person  = p,
        Level   = row_means[p],
        R2      = R2_person[p],
        as.vector(rbind(
          Ori  = U[p, ],
          Mean = colMeans(wmat, na.rm = TRUE),
          SE   = apply(wmat, 2, sd, na.rm = TRUE),
          Lo   = apply(wmat, 2, quantile, la, na.rm = TRUE),
          Hi   = apply(wmat, 2, quantile, ua, na.rm = TRUE)
        ))
      )
    })
    result_B <- do.call(rbind, rows)
    rownames(result_B) <- paste0("#", participants)
  }

  list(weights    = result,
       weights_B  = result_B,
       R2_k       = R2_k,
       R2_mean    = mean(R2_person))
}


## =============================================================================
## 6.  CSV OUTPUT HELPERS
## =============================================================================

write_fapa_results <- function(bca, prefix) {
  K <- bca$K
  for (k in seq_len(K)) {
    fname <- paste0(prefix, "_CoreProfile", k, ".csv")
    df    <- bca$ci[[k]]
    rownames(df) <- bca$varnames
    write.csv(df, file = fname)
    message("Wrote: ", fname)
  }
}

write_verification <- function(pa, pr, tc, prefix, K_pa = NULL) {
  ## Stage 1
  f1 <- paste0(prefix, "_Stage1_PA.csv")
  write.csv(data.frame(Component    = seq_along(pa$obs_sv2),
                       Obs_sv2      = pa$obs_sv2,
                       Obs_PropVar  = pa$prop_obs,
                       Rand_95pct   = pa$thresh,
                       Rand_PropVar = pa$prop_rand,
                       Retained     = pa$obs_sv2 > pa$thresh),
            file = f1, row.names = FALSE)

  ## Stage 2
  f2 <- paste0(prefix, "_Stage2_Procrustes.csv")
  df2 <- data.frame(Dimension   = seq_len(pr$K),
                    Mean_deg    = pr$angle_mean,
                    SD_deg      = pr$angle_sd,
                    CI_lo       = pr$angle_q025,
                    CI_hi       = pr$angle_q975,
                    Angle_OK    = pr$angle_mean < pr$angle_thresh,
                    Prop_stable = pr$prop_stable)
  if (!is.null(K_pa))
    df2$PA_Retained <- seq_len(pr$K) <= K_pa
  write.csv(df2, file = f2, row.names = FALSE)

  ## Stage 3
  f3 <- paste0(prefix, "_Stage3_TuckerCC.csv")
  df3 <- data.frame(Core_Profile = seq_len(tc$K),
                    Mean_CC      = tc$cc_mean,
                    SD_CC        = tc$cc_sd,
                    CI_lo        = tc$cc_q025,
                    CI_hi        = tc$cc_q975,
                    CC_OK        = tc$cc_mean >= tc$cc_thresh)
  if (!is.null(K_pa))
    df3$PA_Retained <- seq_len(tc$K) <= K_pa
  write.csv(df3, file = f3, row.names = FALSE)

  message("Verification results written to: ", f1, ", ", f2, ", ", f3)
}


## =============================================================================
## 7.  PLOTTING FUNCTIONS
## =============================================================================

## Core-profile plot with BCa CIs (split at subscale boundary)
plot_fapa_core <- function(bca, i = 1, split_at = 11,
                            main = NULL, ylab = "Core-Profile Coordinate") {
  if (is.null(main))
    main <- sprintf("FAPA Core Profile %d  (95%% BCa CI)", i)
  s  <- bca$ci[[i]]
  op <- par(mar = c(2.1, 4.1, 2.5, 1)); on.exit(par(op))
  yr <- range(s$BCaLower, s$BCaUpper, na.rm = TRUE)

  idx_before <- seq_len(split_at)
  idx_after  <- seq(split_at + 1L, nrow(s))

  plot(s$Ori[idx_before], type = "b", lty = 2, lwd = 2,
       col = "red",  ylim = yr, xlab = "", ylab = ylab, main = main,
       xaxt = "n")
  axis(1, at = seq_along(idx_before), labels = bca$varnames[idx_before],
       las = 2, cex.axis = 0.6)
  lines(s$BCaLower[idx_before], lty = 3, col = "red")
  lines(s$BCaUpper[idx_before], lty = 3, col = "red")
  abline(h = 0, lty = 3, col = "grey50")
  lines(s$Ori[idx_after], type = "b", lwd = 2, col = "blue")
  lines(s$BCaLower[idx_after], lty = 3, col = "blue")
  lines(s$BCaUpper[idx_after], lty = 3, col = "blue")
  legend("bottom", inset = -0.02, xpd = NA, bty = "n",
         legend = c("Before", "After"),
         col = c("red", "blue"), lty = c(2, 1), lwd = 2, ncol = 2)
}

## PA scree plot: observed σ² vs random 95th-percentile reference line
plot_pa_scree <- function(pa, main = "Horn's Parallel Analysis — Scree") {
  K   <- length(pa$obs_sv2)
  op  <- par(mar = c(4, 4.5, 2.5, 1)); on.exit(par(op))
  ylim <- range(c(pa$obs_sv2, pa$thresh)) * c(0.95, 1.05)
  plot(seq_len(K), pa$obs_sv2, type = "b", lwd = 2, col = "steelblue",
       pch = 16, ylim = ylim,
       main = main, xlab = "Component", ylab = expression(sigma^2))
  lines(seq_len(K), pa$thresh, type = "b", lwd = 2, col = "red",
        lty = 2, pch = 1)
  abline(v = pa$n_retain + 0.5, lty = 3, col = "grey40")
  legend("topright", bty = "n",
         legend = c("Observed", "Random 95th pct", "Retention cut"),
         col    = c("steelblue", "red", "grey40"),
         lty    = c(1, 2, 3), lwd = 2, pch = c(16, 1, NA))
}

## Person–core overlay (standardized)
plot_person_match <- function(bca, Xtilde, p = 1, K = 2) {
  op <- par(mfrow = c(1L, K), mar = c(4, 4.1, 2.5, 1)); on.exit(par(op))
  y  <- scale(as.numeric(Xtilde[p, ]))
  for (i in seq_len(K)) {
    cp <- scale(bca$X0[, i])
    plot(y, type = "b", col = "steelblue",
         main  = sprintf("FAPA CP_%d vs Person #%d", i, p),
         xlab  = "Subscale", ylab = "Standardized")
    lines(cp, lty = 2, col = "tomato", lwd = 2)
    abline(h = 0, lty = 3, col = "grey70")
    legend("topright", bty = "n",
           legend = c("Person", sprintf("CP%d", i)),
           col = c("steelblue","tomato"), lty = c(1,2), lwd = 2)
  }
}

## Principal-angle distribution plot (Stage 2)
plot_principal_angles <- function(pr) {
  K  <- pr$K
  op <- par(mfrow = c(1L, K), mar = c(4, 4.1, 2.5, 1)); on.exit(par(op))
  for (k in seq_len(K)) {
    hist(pr$angles_mat[, k], breaks = 40,
         main  = sprintf("Principal Angle — Dimension %d", k),
         xlab  = "Degrees", col = "lightsteelblue", border = "white")
    abline(v = pr$angle_thresh, col = "red", lty = 2, lwd = 2)
    legend("topright", bty = "n",
           legend = sprintf("< %.0f° stability bound", pr$angle_thresh),
           col = "red", lty = 2, lwd = 2)
  }
}

## Tucker CC distribution plot (Stage 3)
plot_tucker_cc <- function(tc, cc_thresh = 0.85) {
  K  <- tc$K
  op <- par(mfrow = c(1L, K), mar = c(4, 4.1, 2.5, 1)); on.exit(par(op))
  for (k in seq_len(K)) {
    hist(tc$cc_mat[, k], breaks = 40, xlim = c(-1, 1),
         main  = sprintf("Tucker CC — Core Profile %d", k),
         xlab  = "Congruence Coefficient",
         col   = "lightsteelblue", border = "white")
    abline(v = cc_thresh, col = "red",    lty = 2, lwd = 2)
    abline(v = 0.95,      col = "orange", lty = 3, lwd = 2)
    legend("topleft", bty = "n",
           legend = c(sprintf("%.2f (fair)", cc_thresh), "0.95 (high)"),
           col = c("red","orange"), lty = c(2, 3), lwd = 2)
  }
}


## =============================================================================
## 8.  MAIN DRIVER
## =============================================================================

set.seed(seed)
testname <- c(before_labels, after_labels)

## -- Load and ipsatize --------------------------------------------------------
message("Loading data and ipsatizing...")
dat   <- load_and_ipsatize(data_path, testname)
raw   <- dat$raw
Xtilde <- dat$ipsatized
message(sprintf("  %d persons × %d variables", nrow(Xtilde), ncol(Xtilde)))


## -- Stage 1: Horn's Parallel Analysis ----------------------------------------
message("\nRunning Stage 1 — Horn's Parallel Analysis (B = ", B_boot, ")...")
pa_result <- fapa_pa(Xtilde, B = B_boot, alpha = alpha, seed = seed)
print_pa(pa_result)
K_pa     <- pa_result$n_retain                       # dimensionality from PA
K_max    <- length(pa_result$obs_sv2)                # upper bound (I-1 or P-1)
K_report <- min(K_pa + K_extra, K_max)               # extended K for verification display
message(sprintf("  Reporting Stages 2 & 3 for K = %d (retained: %d, extra: %d)",
                K_report, K_pa, K_report - K_pa))


## -- Initial FAPA solution at K_pa components ---------------------------------
message(sprintf("\nFitting FAPA with K = %d component(s)...", K_pa))
fit_fapa <- fapa_core(Xtilde, K = K_pa, direction = direction)

cat(sprintf("\n  Total ipsatized variance     : %.4f\n", fit_fapa$total_var))
cat("  Variance per component      :",
    paste(round(fit_fapa$var_k, 4), collapse = "  "), "\n")
cat("  Proportion explained        :",
    paste(round(fit_fapa$prop_var, 4), collapse = "  "), "\n")
cat("  Cumulative proportion       :",
    paste(round(fit_fapa$cum_var, 4), collapse = "  "), "\n")


## -- Stage 2: Procrustes / Principal Angles -----------------------------------
message("\nRunning Stage 2 — Procrustes / Principal Angles (B = ", B_boot, ")...")
pr_result <- fapa_procrustes(Xtilde, K = K_report, B = B_boot,
                              angle_thresh = angle_thresh, seed = seed)
print_procrustes(pr_result, K_pa)


## -- Stage 3: Tucker's Congruence Coefficients --------------------------------
message("\nRunning Stage 3 — Tucker's Congruence Coefficients (B = ", B_boot, ")...")
tc_result <- fapa_tucker(Xtilde, K = K_report, B = B_boot,
                          cc_thresh = cc_thresh, seed = seed)
print_tucker(tc_result, cc_thresh, K_pa)


## -- BCa confidence intervals for core profiles --------------------------------
message("\nComputing BCa CIs for core profiles (B = ", B_boot, ")...")
bca_result <- fapa_bca(Xtilde, K = K_pa, B = B_boot,
                        alpha = alpha, seed = seed)


## -- Person-level reconstruction ----------------------------------------------
message("\nComputing person-level weights and reconstruction R²...")
person_result <- fapa_person(Xtilde, fit_fapa,
                              participants = participants,
                              B_boot = B_boot, alpha = alpha, seed = seed)

cat(sprintf("\n  Mean person reconstruction R² : %.4f\n",
            person_result$R2_mean))


## -- Correlation of CP1 with subscale grand means (sanity check) --------------
subscale_mean  <- colMeans(raw)
cp1            <- fit_fapa$X[, 1]
cor_cp1_means  <- cor(subscale_mean, cp1)
message(sprintf("\n  Sanity check — cor(grand means, CP1) = %.3f", cor_cp1_means))
if (abs(cor_cp1_means) > 0.70)
  message("  NOTE: CP1 may be confounded with profile level.")


## -- Write outputs ------------------------------------------------------------
message("\nWriting output files...")
write_fapa_results(bca_result, prefix = out_prefix)
write_verification(pa_result, pr_result, tc_result, prefix = out_prefix, K_pa = K_pa)

weights_file <- paste0(out_prefix, "_PersonWeights.csv")
write.csv(person_result$weights, file = weights_file, row.names = FALSE)
message("Wrote: ", weights_file)

if (!is.null(person_result$weights_B)) {
  weights_B_file <- paste0(out_prefix, "_PersonWeights_Bootstrap.csv")
  write.csv(person_result$weights_B, file = weights_B_file)
  message("Wrote: ", weights_B_file)
}


## -- Plots --------------------------------------------------------------------
message("\nGenerating plots...")

## Stage 1: PA scree
plot_pa_scree(pa_result)

## Core profiles with BCa CIs
for (k in seq_len(K_pa))
  plot_fapa_core(bca_result, i = k, split_at = 11)

## Person overlay (first two cores, first participant)
if (K_pa >= 2)
  plot_person_match(bca_result, Xtilde, p = 1, K = min(2, K_pa))

## Verification plots
plot_principal_angles(pr_result)
plot_tucker_cc(tc_result, cc_thresh = cc_thresh)

message("\nFAPA analysis complete.")

## =============================================================================
## END OF fapa_main.R
## =============================================================================
