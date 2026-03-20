## tests/testthat/test-fapa_person.R
## Sanity checks for fapa_person()

test_that("fapa_person returns correct number of rows", {
  set.seed(1)
  P <- 50;  I <- 8;  K <- 2
  Xtilde <- matrix(rnorm(P * I), P, I)
  Xtilde <- Xtilde - rowMeans(Xtilde)
  fit    <- fapa_core(Xtilde, K)

  res <- fapa_person(Xtilde, fit, participants = NULL)
  expect_equal(nrow(res$weights), P)
})

test_that("fapa_person R2 is in [0, 1]", {
  set.seed(2)
  P <- 60;  I <- 10;  K <- 3
  Xtilde <- matrix(rnorm(P * I), P, I)
  Xtilde <- Xtilde - rowMeans(Xtilde)
  fit    <- fapa_core(Xtilde, K)

  res <- fapa_person(Xtilde, fit)
  expect_true(all(res$weights$R2 >= -1e-8))
  expect_true(all(res$weights$R2 <=  1 + 1e-8))
})

test_that("fapa_person weight columns match K", {
  set.seed(3)
  P <- 40;  I <- 8;  K <- 2
  Xtilde <- matrix(rnorm(P * I), P, I)
  Xtilde <- Xtilde - rowMeans(Xtilde)
  fit    <- fapa_core(Xtilde, K)

  res   <- fapa_person(Xtilde, fit)
  wcols <- paste0("w", seq_len(K))
  expect_true(all(wcols %in% names(res$weights)))
})

test_that("fapa_person bootstrap output has correct participant rows", {
  set.seed(4)
  P    <- 50;  I <- 8;  K <- 2
  Xtilde <- matrix(rnorm(P * I), P, I)
  Xtilde <- Xtilde - rowMeans(Xtilde)
  fit    <- fapa_core(Xtilde, K)
  parts  <- c(1L, 5L, 10L)

  res <- fapa_person(Xtilde, fit, participants = parts,
                     B_boot = 30, seed = 4)
  expect_equal(nrow(res$weights_B), length(parts))
})
