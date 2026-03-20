## tests/testthat/test-fapa_core.R
## Basic sanity checks for fapa_core() and load_and_ipsatize()

test_that("fapa_core returns correct dimensions", {
  set.seed(1)
  P <- 50;  I <- 10;  K <- 2
  Xtilde <- matrix(rnorm(P * I), P, I)
  Xtilde <- Xtilde - rowMeans(Xtilde)   # ipsatize

  fit <- fapa_core(Xtilde, K)

  expect_equal(dim(fit$U), c(P, K))
  expect_equal(dim(fit$V), c(I, K))
  expect_equal(dim(fit$X), c(I, K))
  expect_length(fit$S, K)
})

test_that("fapa_core exact reconstruction holds", {
  set.seed(2)
  P <- 40;  I <- 8;  K <- 3
  Xtilde <- matrix(rnorm(P * I), P, I)
  Xtilde <- Xtilde - rowMeans(Xtilde)

  fit  <- fapa_core(Xtilde, K)
  Xhat <- fit$U %*% t(fit$X)
  res  <- Xtilde - Xhat

  ## Residual should equal sum of remaining singular components
  sv_full <- svd(Xtilde, nu = 0, nv = 0)
  expected_rss <- sum(sv_full$d[(K + 1):min(P - 1, I - 1)]^2)
  expect_equal(sum(res^2), expected_rss, tolerance = 1e-6)
})

test_that("fapa_core proportions sum to <= 1", {
  set.seed(3)
  P <- 60;  I <- 12;  K <- 4
  Xtilde <- matrix(rnorm(P * I), P, I)
  Xtilde <- Xtilde - rowMeans(Xtilde)

  fit <- fapa_core(Xtilde, K)
  expect_lte(fit$cum_var[K], 1 + 1e-10)
  expect_true(all(fit$prop_var >= 0))
})

test_that("load_and_ipsatize row sums are zero", {
  tmp <- tempfile(fileext = ".csv")
  dat <- as.data.frame(matrix(sample(1:5, 30, replace = TRUE), nrow = 6))
  write.csv(dat, tmp, row.names = FALSE)

  res <- load_and_ipsatize(tmp, col_labels = paste0("V", 1:5))
  row_sums <- rowSums(res$ipsatized)
  expect_equal(max(abs(row_sums)), 0, tolerance = 1e-10)

  unlink(tmp)
})
