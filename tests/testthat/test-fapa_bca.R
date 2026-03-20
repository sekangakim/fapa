## tests/testthat/test-fapa_bca.R
## Sanity checks for fapa_bca()

test_that("fapa_bca returns K CI data frames with correct columns", {
  set.seed(1)
  P <- 60;  I <- 8;  K <- 2
  Xtilde <- matrix(rnorm(P * I), P, I)
  Xtilde <- Xtilde - rowMeans(Xtilde)

  bca <- fapa_bca(Xtilde, K = K, B = 50, seed = 1)

  expect_length(bca$ci, K)
  expected_cols <- c("Ori", "Mean", "SE", "Lower", "Upper",
                     "BCaLower", "BCaUpper")
  for (k in seq_len(K))
    expect_named(bca$ci[[k]], expected_cols)
})

test_that("fapa_bca CI bounds are ordered: BCaLower <= Ori <= BCaUpper", {
  set.seed(2)
  P <- 80;  I <- 10;  K <- 2
  Xtilde <- matrix(rnorm(P * I), P, I)
  Xtilde <- Xtilde - rowMeans(Xtilde)

  bca <- fapa_bca(Xtilde, K = K, B = 100, seed = 2)

  for (k in seq_len(K)) {
    df <- bca$ci[[k]]
    ## BCaLower should be <= Ori (allowing small numerical noise)
    expect_true(all(df$BCaLower <= df$Ori + 1e-6))
    ## BCaUpper should be >= Ori
    expect_true(all(df$BCaUpper >= df$Ori - 1e-6))
  }
})

test_that("fapa_bca boot_X array has correct dimensions", {
  set.seed(3)
  P <- 50;  I <- 8;  K <- 2;  B <- 40
  Xtilde <- matrix(rnorm(P * I), P, I)
  Xtilde <- Xtilde - rowMeans(Xtilde)

  bca <- fapa_bca(Xtilde, K = K, B = B, seed = 3)
  expect_equal(dim(bca$boot_X), c(B, I, K))
})

test_that(".align_signs corrects sign flips", {
  set.seed(4)
  I <- 10;  K <- 2
  X0 <- matrix(rnorm(I * K), I, K)
  Xb <- X0 * c(1, -1)               # flip second column

  Xb_aligned <- FAPA:::.align_signs(X0, Xb)
  expect_equal(Xb_aligned, X0, tolerance = 1e-12)
})
