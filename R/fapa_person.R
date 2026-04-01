# =============================================================================
# fapa_person.R
# Person-level reconstruction, pattern weights, and bootstrap CIs for
# selected participants.
# =============================================================================


#' Person-level reconstruction and pattern weights
#'
#' Projects each person onto the retained core profiles, returning
#' reconstruction \eqn{R^2}, pattern weights \eqn{\mathbf{u}_p}, and
#' normalised person-core correlations.  Optionally computes percentile
#' bootstrap confidence intervals for pattern weights of selected participants.
#'
#' @param Xtilde       Numeric matrix (persons \eqn{\times} variables), already
#'   ipsatized.
#' @param fit          A list returned by \code{\link{fapa_core}}.
#' @param participants Integer vector of row indices for which bootstrap CIs
#'   are desired, or \code{NULL} to skip individual inference.
#' @param B_boot       Integer. Bootstrap replicates for participant CIs.
#'   Default \code{2000}.
#' @param alpha        Numeric. Two-tailed significance level. Default
#'   \code{0.05}.
#' @param seed         Integer or \code{NULL}. Random seed.
#'
#' @return A named list:
#'   \describe{
#'     \item{weights}{Data frame (\eqn{P} rows) with columns \code{Person},
#'       \code{Level}, \code{R2}, \code{w1} \ldots \code{wK},
#'       \code{rDim1} \ldots \code{rDimK}.}
#'     \item{weights_B}{Matrix of bootstrap summary statistics for
#'       \code{participants} (or \code{NULL}).}
#'     \item{R2_k}{Partial \eqn{R^2} per dimension.}
#'     \item{R2_mean}{Mean person reconstruction \eqn{R^2} across all persons.}
#'   }
#'
#' @export
fapa_person <- function(Xtilde, fit, participants = NULL,
                         B_boot = 2000, alpha = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  P  <- nrow(Xtilde)
  K  <- fit$K
  la <- alpha / 2
  ua <- 1 - alpha / 2

  U         <- fit$U
  X         <- fit$X
  row_means <- rowMeans(Xtilde)

  ## Partial R-squared per dimension
  R2_k <- numeric(K)
  for (k in seq_len(K))
    R2_k[k] <- 1 - sum((Xtilde - U[, k, drop = FALSE] %*%
                           t(X[, k, drop = FALSE]))^2) / sum(Xtilde^2)

  ## Full reconstruction R-squared per person
  Xhat      <- U %*% t(X)
  R2_person <- 1 - rowSums((Xtilde - Xhat)^2) / rowSums(Xtilde^2)

  result <- data.frame(
    Person = seq_len(P),
    Level  = row_means,
    R2     = round(R2_person, 4),
    matrix(round(U,              4), ncol = K,
           dimnames = list(NULL, paste0("w",    seq_len(K)))),
    matrix(round(fit$person_cor, 4), ncol = K,
           dimnames = list(NULL, paste0("rDim", seq_len(K))))
  )

  ## Bootstrap CIs for selected participants
  result_B <- NULL
  if (!is.null(participants)) {
    boot_w <- array(NA_real_, dim = c(B_boot, length(participants), K))

    for (b in seq_len(B_boot)) {
      idx  <- sample(P, P, replace = TRUE)
      Xb   <- Xtilde[idx, , drop = FALSE]
      fitb <- fapa_core(Xb, K, direction = fit$direction)
      for (pi in seq_along(participants)) {
        yp <- as.numeric(Xtilde[participants[pi], ])
        boot_w[b, pi, ] <- as.numeric(
          solve(crossprod(fitb$X), crossprod(fitb$X, yp)))
      }
    }

    rows <- lapply(seq_along(participants), function(pi) {
      p    <- participants[pi]
      wmat <- boot_w[, pi, , drop = FALSE]
      wmat <- matrix(wmat, nrow = B_boot, ncol = K)
      c(Person = p,
        Level  = row_means[p],
        R2     = R2_person[p],
        as.vector(rbind(
          Ori  = U[p, ],
          Mean = colMeans(wmat,                           na.rm = TRUE),
          SE   = apply(wmat, 2, stats::sd,                na.rm = TRUE),
          Lo   = apply(wmat, 2, stats::quantile, la,      na.rm = TRUE),
          Hi   = apply(wmat, 2, stats::quantile, ua,      na.rm = TRUE)
        )))
    })
    result_B <- do.call(rbind, rows)
    rownames(result_B) <- paste0("#", participants)
  }

  list(weights   = result,
       weights_B = result_B,
       R2_k      = R2_k,
       R2_mean   = mean(R2_person))
}
