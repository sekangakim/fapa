# =============================================================================
# fapa_bca.R
# BCa bootstrap confidence intervals for core-profile coordinates.
# Uses boot::boot() + boot::boot.ci(type = "bca") — Davison & Hinkley (1997).
# =============================================================================


#' Align bootstrap core-profile columns to a reference via inner-product signs
#'
#' For each column \eqn{k}, flips the sign of \code{Xb[, k]} if its inner
#' product with \code{X0[, k]} is negative.  This prevents the bootstrap
#' distribution of core-profile coordinates from becoming bimodal due to SVD
#' sign ambiguity across resamples.
#'
#' @param X0 Reference core-profile matrix (\eqn{I \times K}).
#' @param Xb Bootstrap core-profile matrix (\eqn{I \times K}).
#' @return \code{Xb} with columns sign-aligned to \code{X0}.
#' @keywords internal
.align_signs <- function(X0, Xb) {
  for (k in seq_len(ncol(X0)))
    if (sum(X0[, k] * Xb[, k]) < 0) Xb[, k] <- -Xb[, k]
  Xb
}


#' BCa bootstrap confidence intervals for FAPA core profiles
#'
#' Computes BCa (bias-corrected and accelerated) bootstrap confidence intervals
#' for every coordinate of every retained core profile, using the canonical
#' implementation in \code{\link[boot]{boot}} and \code{\link[boot]{boot.ci}}.
#'
#' Sign ambiguity across bootstrap resamples is handled inside the statistic
#' function via an inner-product alignment rule (see \code{.align_signs}),
#' ensuring that each bootstrap distribution is unimodal before BCa adjustment.
#'
#' @param Xtilde Numeric matrix (persons × variables), already ipsatized.
#' @param K      Integer. Number of core profiles (must equal the retained
#'   dimensionality from \code{\link{fapa_pa}}).
#' @param B      Integer. Number of bootstrap replicates. Default \code{2000}.
#' @param alpha  Numeric. Two-tailed significance level. Default \code{0.05}.
#' @param seed   Integer or \code{NULL}. Random seed.
#'
#' @return A named list:
#'   \describe{
#'     \item{ci}{List of \eqn{K} data frames (one per core profile), each with
#'       columns \code{Ori}, \code{Mean}, \code{SE}, \code{Lower},
#'       \code{Upper}, \code{BCaLower}, \code{BCaUpper}.}
#'     \item{X0}{Original core-profile matrix (\eqn{I \times K}).}
#'     \item{boot_out}{The full \code{boot} object for downstream diagnostics.}
#'     \item{boot_X}{3-D array (\eqn{B \times I \times K}) of bootstrap profiles.}
#'     \item{K, B, alpha, varnames}{Inputs echoed for plotting and output.}
#'   }
#'
#' @references
#' Davison, A. C., & Hinkley, D. V. (1997).
#' \emph{Bootstrap Methods and Their Application}.
#' Cambridge University Press. \doi{10.1017/CBO9780511802843}
#'
#' @seealso \code{\link{plot_fapa_core}}, \code{\link{write_fapa_results}}
#'
#' @export
fapa_bca <- function(Xtilde, K, B = 2000, alpha = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  I <- ncol(Xtilde)

  fit0 <- fapa_core(Xtilde, K)
  X0   <- fit0$X                              # I × K original core profiles

  ## Statistic: returns aligned core profiles as a flat vector (length I*K)
  fapa_stat <- function(data, indices) {
    Xb <- fapa_core(data[indices, , drop = FALSE], K)$X
    as.vector(.align_signs(X0, Xb))
  }

  boot_out <- boot::boot(data = Xtilde, statistic = fapa_stat, R = B)

  la      <- alpha / 2
  ua      <- 1 - alpha / 2
  n_stats <- I * K

  bca_lo  <- numeric(n_stats);  bca_hi  <- numeric(n_stats)
  boot_mn <- numeric(n_stats);  boot_se <- numeric(n_stats)
  pct_lo  <- numeric(n_stats);  pct_hi  <- numeric(n_stats)

  for (idx in seq_len(n_stats)) {
    ci <- tryCatch(
      boot::boot.ci(boot_out, conf = 1 - alpha, type = "bca", index = idx),
      error = function(e) NULL
    )
    if (!is.null(ci)) {
      bca_lo[idx] <- ci$bca[4]
      bca_hi[idx] <- ci$bca[5]
    } else {
      ## Fallback to percentile when BCa cannot be computed
      bca_lo[idx] <- quantile(boot_out$t[, idx], la, na.rm = TRUE)
      bca_hi[idx] <- quantile(boot_out$t[, idx], ua, na.rm = TRUE)
    }
    boot_mn[idx] <- mean(boot_out$t[, idx],     na.rm = TRUE)
    boot_se[idx] <- sd(boot_out$t[, idx],        na.rm = TRUE)
    pct_lo[idx]  <- quantile(boot_out$t[, idx], la, na.rm = TRUE)
    pct_hi[idx]  <- quantile(boot_out$t[, idx], ua, na.rm = TRUE)
  }

  ci_list <- vector("list", K)
  for (k in seq_len(K)) {
    idx_k <- seq_len(I) + (k - 1L) * I
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

  list(ci       = ci_list,
       X0       = X0,
       boot_out = boot_out,
       boot_X   = array(boot_out$t, dim = c(B, I, K)),
       K        = K,
       B        = B,
       alpha    = alpha,
       varnames = colnames(Xtilde))
}
