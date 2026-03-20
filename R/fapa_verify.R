# =============================================================================
# fapa_verify.R
# Three-stage bootstrap verification framework:
#   Stage 1 – Horn's parallel analysis        (dimensionality)
#   Stage 2 – Procrustes / principal angles   (subspace stability)
#   Stage 3 – Tucker's congruence coefficient (profile replicability)
# =============================================================================


# -----------------------------------------------------------------------------
# Internal helpers (not exported)
# -----------------------------------------------------------------------------

#' Compute principal angles between two column spaces (degrees)
#'
#' @param A Numeric matrix whose columns span the first subspace.
#' @param B Numeric matrix whose columns span the second subspace.
#'   Both matrices must have the same number of columns; they need not be
#'   orthonormal.
#' @return Numeric vector of principal angles in degrees (length = ncol(A)).
#' @keywords internal
principal_angles_deg <- function(A, B) {
  QA      <- qr.Q(qr(A))
  QB      <- qr.Q(qr(B))
  cos_ang <- pmin(svd(crossprod(QA, QB), nu = 0, nv = 0)$d, 1)
  acos(cos_ang) * 180 / pi
}


#' Tucker's congruence coefficient between two vectors
#'
#' Computes the cosine of the angle between \code{x} and \code{y}, equivalent
#' to Tucker's (1951) congruence coefficient.
#'
#' @param x,y Numeric vectors of the same length.
#' @return Scalar; \code{NA} if either vector has zero norm.
#' @keywords internal
tucker_cc <- function(x, y) {
  denom <- sqrt(sum(x^2) * sum(y^2))
  if (denom < .Machine$double.eps) return(NA_real_)
  sum(x * y) / denom
}


# -----------------------------------------------------------------------------
# Stage 1 – Horn's Parallel Analysis
# -----------------------------------------------------------------------------

#' Stage 1: Horn's Parallel Analysis (variance-matched bootstrap)
#'
#' Determines the number of components to retain from the SVD of an ipsatized
#' data matrix using a variance-matched bootstrap version of Horn's (1965)
#' parallel analysis.
#'
#' For each of \code{B} bootstrap replicates, a random matrix of identical
#' dimensions is row-centred (ipsatized) and then rescaled to the same
#' Frobenius norm as \code{Xtilde}.  This variance-matching step is essential:
#' without it, raw-score data trivially dominates N(0,1) random matrices and
#' PA retains all components.  Components whose observed \eqn{\sigma^2_k}
#' exceeds the \eqn{(1-\alpha)} quantile of the matched null distribution are
#' retained.
#'
#' @param Xtilde Numeric matrix (persons × variables), already ipsatized.
#' @param B      Integer. Number of bootstrap replicates. Default \code{2000}.
#' @param alpha  Numeric. Significance level. Default \code{0.05}.
#' @param seed   Integer or \code{NULL}. Random seed for reproducibility.
#'
#' @return A named list:
#'   \describe{
#'     \item{n_retain}{Number of components retained.}
#'     \item{obs_sv2}{Observed squared singular values (length = \eqn{K_{\max}}).}
#'     \item{thresh}{Bootstrap \eqn{(1-\alpha)} quantile per component.}
#'     \item{prop_obs}{Proportion of variance per observed component.}
#'     \item{prop_rand}{Mean proportion of variance per random component.}
#'     \item{rand_sv2}{Full \eqn{B \times K_{\max}} matrix of random \eqn{\sigma^2}.}
#'     \item{total_var}{Total ipsatized variance.}
#'     \item{alpha, B}{Inputs echoed for reporting.}
#'   }
#'
#' @references
#' Horn, J. L. (1965). A rationale and test for the number of factors in
#' factor analysis. \emph{Psychometrika, 30}(2), 179--185.
#'
#' @seealso \code{\link{print_pa}}, \code{\link{plot_pa_scree}}
#'
#' @export
fapa_pa <- function(Xtilde, B = 2000, alpha = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  P     <- nrow(Xtilde);  I <- ncol(Xtilde)
  K_max <- min(P - 1L, I - 1L)

  obs_sv2   <- svd(Xtilde, nu = 0, nv = 0)$d[seq_len(K_max)]^2
  total_var <- sum(obs_sv2)

  rand_sv2 <- matrix(0, B, K_max)
  for (b in seq_len(B)) {
    Xr <- matrix(rnorm(P * I), P, I)
    Xr <- Xr - rowMeans(Xr)                           # ipsatize
    Xr <- Xr * sqrt(total_var / sum(Xr^2))            # variance-match
    rand_sv2[b, ] <- svd(Xr, nu = 0, nv = 0)$d[seq_len(K_max)]^2
  }

  thresh   <- apply(rand_sv2, 2, quantile, probs = 1 - alpha)
  n_retain <- sum(obs_sv2 > thresh)
  if (n_retain == 0L) {
    warning("PA retained 0 components. Check data scaling or reduce alpha.")
    n_retain <- 1L
  }

  list(n_retain  = n_retain,
       obs_sv2   = obs_sv2,
       thresh    = thresh,
       prop_obs  = obs_sv2 / total_var,
       prop_rand = colMeans(rand_sv2) / total_var,
       rand_sv2  = rand_sv2,
       total_var = total_var,
       alpha     = alpha,
       B         = B)
}


#' Print a summary of Stage 1 parallel analysis results
#'
#' @param pa A list returned by \code{\link{fapa_pa}}.
#' @return Invisibly returns \code{NULL}. Called for its side-effect of
#'   printing a formatted table to the console.
#' @export
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
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Stage 2 – Procrustes / Principal Angles
# -----------------------------------------------------------------------------

#' Stage 2: Procrustes stability via principal angles
#'
#' For each of \code{B} bootstrap resamples of the ipsatized data matrix,
#' computes the \eqn{K} right singular vectors and measures the principal
#' angles (in degrees) between the bootstrap subspace and the original
#' \eqn{K}-dimensional right singular vector subspace.
#'
#' A bootstrap replicate is declared \emph{stable} when \strong{all}
#' \eqn{K} principal angles are strictly less than \code{angle_thresh}.
#' This criterion confirms that the bootstrap subspace is nearly parallel to
#' the original, providing geometric evidence of dimensional stability.
#'
#' @param Xtilde       Numeric matrix (persons × variables), already ipsatized.
#' @param K            Integer. Number of dimensions to assess.
#' @param B            Integer. Number of bootstrap replicates. Default \code{2000}.
#' @param angle_thresh Numeric. Upper stability bound in degrees. Default \code{30}.
#' @param seed         Integer or \code{NULL}. Random seed.
#'
#' @return A named list:
#'   \describe{
#'     \item{angles_mat}{B × K matrix of principal angles (degrees).}
#'     \item{angle_mean, angle_sd}{Per-dimension mean and SD of angles.}
#'     \item{angle_q025, angle_q975}{Per-dimension 2.5th and 97.5th percentiles.}
#'     \item{n_stable}{Number of replicates satisfying the stability criterion.}
#'     \item{prop_stable}{Proportion of stable replicates.}
#'     \item{angle_thresh, K, B}{Inputs echoed for reporting.}
#'   }
#'
#' @references
#' Björck, Å., & Golub, G. H. (1973). Numerical methods for computing angles
#' between linear subspaces. \emph{Mathematics of Computation, 27}(123),
#' 579--594.
#'
#' @seealso \code{\link{print_procrustes}}, \code{\link{plot_principal_angles}}
#'
#' @export
fapa_procrustes <- function(Xtilde, K, B = 2000,
                             angle_thresh = 30, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n  <- nrow(Xtilde)
  V0 <- svd(Xtilde, nu = 0, nv = K)$v          # I × K original subspace

  angles_mat <- matrix(NA_real_, B, K)
  n_stable   <- 0L

  for (b in seq_len(B)) {
    idx <- sample(n, n, replace = TRUE)
    Vb  <- svd(Xtilde[idx, , drop = FALSE], nu = 0, nv = K)$v
    ang <- principal_angles_deg(V0, Vb)
    angles_mat[b, ] <- ang
    if (all(ang < angle_thresh)) n_stable <- n_stable + 1L
  }

  list(angles_mat   = angles_mat,
       angle_mean   = colMeans(angles_mat),
       angle_sd     = apply(angles_mat, 2, sd),
       angle_q025   = apply(angles_mat, 2, quantile, 0.025),
       angle_q975   = apply(angles_mat, 2, quantile, 0.975),
       n_stable     = n_stable,
       prop_stable  = n_stable / B,
       angle_thresh = angle_thresh,
       K = K, B = B)
}


#' Print a summary of Stage 2 principal-angle results
#'
#' @param pr   A list returned by \code{\link{fapa_procrustes}}.
#' @param K_pa Integer or \code{NULL}. If supplied, a \code{PA_Retained}
#'   column is appended showing which dimensions were retained by PA.
#' @return Invisibly returns \code{NULL}.
#' @export
print_procrustes <- function(pr, K_pa = NULL) {
  cat("\n=== STAGE 2: Procrustes / Principal Angles ===\n")
  cat(sprintf("  Bootstrap replicates           : %d\n", pr$B))
  cat(sprintf("  Stability threshold            : < %.0f degrees\n",
              pr$angle_thresh))
  cat(sprintf(
    "  Replicates with ALL angles < threshold: %d (%.1f%%)\n\n",
    pr$n_stable, 100 * pr$prop_stable))
  df <- data.frame(
    Dimension = seq_len(pr$K),
    Mean_deg  = round(pr$angle_mean, 2),
    SD_deg    = round(pr$angle_sd,   2),
    CI_lo_deg = round(pr$angle_q025, 2),
    CI_hi_deg = round(pr$angle_q975, 2),
    Angle_OK  = ifelse(pr$angle_mean < pr$angle_thresh, "YES", "no")
  )
  if (!is.null(K_pa))
    df$PA_Retained <- ifelse(seq_len(pr$K) <= K_pa, "YES", "no")
  print(df, row.names = FALSE)
  invisible(NULL)
}


# -----------------------------------------------------------------------------
# Stage 3 – Tucker's Congruence Coefficients
# -----------------------------------------------------------------------------

#' Stage 3: Profile replicability via Tucker's congruence coefficients
#'
#' For each of \code{B} bootstrap resamples, computes Tucker's congruence
#' coefficient (CC) between each original core profile and its bootstrap
#' counterpart.  Sign ambiguity is resolved by choosing the sign that
#' maximises the absolute CC before storing.
#'
#' Conventional thresholds (Lorenzo-Seva & ten Berge, 2006):
#' \itemize{
#'   \item CC \eqn{\ge} 0.95: high similarity ("factor replication").
#'   \item CC \eqn{\ge} 0.85: fair similarity.
#'   \item CC \eqn{<} 0.85: poor similarity.
#' }
#'
#' @param Xtilde    Numeric matrix (persons × variables), already ipsatized.
#' @param K         Integer. Number of core profiles to assess.
#' @param B         Integer. Number of bootstrap replicates. Default \code{2000}.
#' @param cc_thresh Numeric. Acceptability lower bound. Default \code{0.85}.
#' @param seed      Integer or \code{NULL}. Random seed.
#'
#' @return A named list:
#'   \describe{
#'     \item{cc_mat}{B × K matrix of Tucker CCs.}
#'     \item{cc_mean, cc_sd}{Per-profile mean and SD of CCs.}
#'     \item{cc_q025, cc_q975}{Per-profile 2.5th and 97.5th percentiles.}
#'     \item{cc_thresh, K, B}{Inputs echoed for reporting.}
#'   }
#'
#' @references
#' Lorenzo-Seva, U., & ten Berge, J. M. F. (2006). Tucker's congruence
#' coefficient as a meaningful index of factor similarity.
#' \emph{Methodology, 2}(2), 57--64. \doi{10.1027/1614-2241.2.2.57}
#'
#' Tucker, L. R. (1951). \emph{A method for synthesis of factor analysis
#' studies} (Personnel Research Section Report No. 984). Department of the
#' Army.
#'
#' @seealso \code{\link{print_tucker}}, \code{\link{plot_tucker_cc}}
#'
#' @export
fapa_tucker <- function(Xtilde, K, B = 2000,
                         cc_thresh = 0.85, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n    <- nrow(Xtilde)
  fit0 <- fapa_core(Xtilde, K)
  X0   <- fit0$X                              # I × K original core profiles

  cc_mat <- matrix(NA_real_, B, K)

  for (b in seq_len(B)) {
    idx   <- sample(n, n, replace = TRUE)
    Xb_cp <- fapa_core(Xtilde[idx, , drop = FALSE], K)$X
    for (k in seq_len(K)) {
      cc_pos <- tucker_cc(X0[, k],  Xb_cp[, k])
      cc_neg <- tucker_cc(X0[, k], -Xb_cp[, k])
      cc_mat[b, k] <- ifelse(abs(cc_pos) >= abs(cc_neg), cc_pos, cc_neg)
    }
  }

  list(cc_mat    = cc_mat,
       cc_mean   = colMeans(cc_mat, na.rm = TRUE),
       cc_sd     = apply(cc_mat, 2, sd,       na.rm = TRUE),
       cc_q025   = apply(cc_mat, 2, quantile, 0.025, na.rm = TRUE),
       cc_q975   = apply(cc_mat, 2, quantile, 0.975, na.rm = TRUE),
       cc_thresh = cc_thresh,
       K = K, B = B)
}


#' Print a summary of Stage 3 Tucker CC results
#'
#' @param tc       A list returned by \code{\link{fapa_tucker}}.
#' @param cc_thresh Numeric. Acceptability cutoff to display (should match the
#'   value used in \code{\link{fapa_tucker}}).
#' @param K_pa     Integer or \code{NULL}. If supplied, a \code{PA_Retained}
#'   column is appended.
#' @return Invisibly returns \code{NULL}.
#' @export
print_tucker <- function(tc, cc_thresh, K_pa = NULL) {
  cat("\n=== STAGE 3: Tucker's Congruence Coefficients ===\n")
  cat(sprintf("  Bootstrap replicates : %d\n", tc$B))
  cat(sprintf("  Acceptability cutoff : CC >= %.2f\n\n", cc_thresh))
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
  cat("    CC >= 0.95 = high; CC >= 0.85 = fair; CC < 0.85 = poor\n")
  invisible(NULL)
}
