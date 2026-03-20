## tests/testthat/helper-setup.R
## Shared helpers available to all test files.

#' Generate a small ipsatized matrix for testing
#'
#' @param P Number of persons.
#' @param I Number of variables.
#' @param seed Random seed.
#' @return Ipsatized numeric matrix (P x I).
make_ipsatized <- function(P = 80, I = 10, seed = 42) {
  set.seed(seed)
  X <- matrix(rnorm(P * I), P, I)
  X - rowMeans(X)
}
