# =============================================================================
# fapa_core.R
# Data loading, ipsatization, and core FAPA estimation via SVD.
# =============================================================================


#' Load a CSV and ipsatize (row-centre) it
#'
#' Reads a person-by-variable CSV file, assigns column labels, and removes
#' each person's mean across variables (ipsatization), isolating within-person
#' pattern structure from overall profile elevation.
#'
#' @param path      Character. Path to a CSV file with a header row.
#' @param col_labels Character vector of length equal to the number of columns.
#'   Column names are replaced with these labels after loading.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{raw}{Original data as a \code{data.frame}.}
#'     \item{ipsatized}{Row-centred matrix (\eqn{\tilde{M}}).}
#'     \item{row_means}{Numeric vector of per-person means (profile levels).}
#'     \item{varnames}{The supplied \code{col_labels}.}
#'   }
#'
#' @examples
#' ## Create a small temporary CSV and ipsatize it
#' tmp <- tempfile(fileext = ".csv")
#' write.csv(matrix(sample(1:5, 30, replace = TRUE), nrow = 6),
#'           tmp, row.names = FALSE)
#' dat <- load_and_ipsatize(tmp, col_labels = paste0("V", 1:5))
#' round(rowSums(dat$ipsatized), 10)  # should all be 0
#' unlink(tmp)
#'
#' @export
load_and_ipsatize <- function(path, col_labels) {
  dat <- utils::read.csv(path, header = TRUE, check.names = FALSE)
  if (ncol(dat) != length(col_labels))
    stop(sprintf("Expected %d columns, got %d.", length(col_labels), ncol(dat)))
  colnames(dat) <- col_labels
  row_means <- rowMeans(dat, na.rm = TRUE)
  Xtilde    <- as.matrix(sweep(dat, 1, row_means, FUN = "-"))
  list(raw       = dat,
       ipsatized = Xtilde,
       row_means = row_means,
       varnames  = col_labels)
}


#' Core FAPA estimation via SVD
#'
#' Computes the thin singular value decomposition (SVD) of an ipsatized
#' person-by-variable matrix, returning the \eqn{K} core profiles, person
#' weights, singular values, and variance-accounting summaries.
#'
#' The core-profile (scale) matrix is defined as
#' \eqn{\mathbf{X} = \mathbf{V}\boldsymbol{\Sigma}}, so that each
#' individual's ipsatized profile satisfies
#' \eqn{\tilde{\mathbf{m}}_p = \mathbf{X} \mathbf{u}_p} exactly
#' (rank-\eqn{K} reconstruction).
#'
#' Signs of the singular vectors are normalised so that the element with the
#' largest absolute value in each core-profile column is positive.
#'
#' @param Xtilde   Numeric matrix (persons \eqn{\times} variables), already
#'   ipsatized.
#' @param K        Integer. Number of components to extract.
#' @param direction Integer vector of length \eqn{K} with values \eqn{\pm 1},
#'   overriding the automatic sign convention.  \code{NULL} (default) triggers
#'   automatic normalisation.
#'
#' @return A named list:
#'   \describe{
#'     \item{U}{Person-weight matrix (\eqn{P \times K}).}
#'     \item{S}{Singular values (length \eqn{K}).}
#'     \item{V}{Right singular vectors (\eqn{I \times K}).}
#'     \item{X}{Core-profile matrix \eqn{\mathbf{V}\boldsymbol{\Sigma}}
#'       (\eqn{I \times K}).}
#'     \item{total_var}{Total ipsatized variance (Frobenius norm\eqn{^2}).}
#'     \item{var_k}{Variance per component (\eqn{\sigma_k^2}).}
#'     \item{prop_var}{Proportion of variance per component.}
#'     \item{cum_var}{Cumulative proportion of variance.}
#'     \item{person_cor}{Normalised person-core correlations (\eqn{P \times K}).}
#'     \item{direction}{Sign vector applied (for reproducibility).}
#'     \item{K}{Number of components extracted.}
#'   }
#'
#' @references
#' Kim, S.-K. (2024). Factorization of person response profiles to identify
#' summative profiles carrying central response patterns.
#' \emph{Psychological Methods, 29}(4), 723--730.
#' \doi{10.1037/met0000568}
#'
#' @export
fapa_core <- function(Xtilde, K, direction = NULL) {

  sv  <- svd(Xtilde, nu = K, nv = K)
  U   <- sv$u                            # P x K
  S   <- sv$d[seq_len(K)]               # K singular values
  V   <- sv$v                            # I x K

  ## Core-profile (scale) matrix  X = V Sigma
  X   <- V %*% diag(S, nrow = K)

  ## Sign normalisation: flip so largest-abs element per column is positive
  if (is.null(direction)) {
    direction <- sign(apply(X, 2, function(x) x[which.max(abs(x))]))
    direction[direction == 0] <- 1L
  }
  X <- X %*% diag(direction, nrow = K)
  U <- U %*% diag(direction, nrow = K)
  V <- V %*% diag(direction, nrow = K)

  ## Variance accounting
  total_var <- sum(Xtilde^2)
  var_k     <- S^2
  prop_var  <- var_k / total_var
  cum_var   <- cumsum(prop_var)

  ## Normalised person-core correlations
  prof_norms <- sqrt(rowSums(Xtilde^2))
  person_cor <- sweep(U * rep(S, each = nrow(U)), 1,
                      pmax(prof_norms, .Machine$double.eps), "/")

  list(U          = U,
       S          = S,
       V          = V,
       X          = X,
       total_var  = total_var,
       var_k      = var_k,
       prop_var   = prop_var,
       cum_var    = cum_var,
       person_cor = person_cor,
       direction  = direction,
       K          = K)
}
