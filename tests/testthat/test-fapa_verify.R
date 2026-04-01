## tests/testthat/test-fapa_verify.R
## Sanity checks for fapa_pa(), fapa_procrustes(), fapa_tucker()

## ---- helpers ----------------------------------------------------------------

make_ipsatized <- function(P = 80, I = 10, seed = 42) {
  set.seed(seed)
  X <- matrix(rnorm(P * I), P, I)
  X - rowMeans(X)
}

## ---- Stage 1: fapa_pa -------------------------------------------------------

test_that("fapa_pa retains at least 1 component", {
  Xtilde <- make_ipsatized()
  pa <- fapa_pa(Xtilde, B = 50, seed = 1)
  expect_gte(pa$n_retain, 1L)
})

test_that("fapa_pa output vectors have length K_max", {
  Xtilde <- make_ipsatized(P = 60, I = 8)
  pa     <- fapa_pa(Xtilde, B = 50, seed = 1)
  K_max  <- min(nrow(Xtilde) - 1L, ncol(Xtilde) - 1L)
  expect_length(pa$obs_sv2,   K_max)
  expect_length(pa$thresh,    K_max)
  expect_length(pa$prop_obs,  K_max)
  expect_length(pa$prop_rand, K_max)
})

test_that("fapa_pa proportions are in [0, 1]", {
  Xtilde <- make_ipsatized()
  pa     <- fapa_pa(Xtilde, B = 50, seed = 1)
  expect_true(all(pa$prop_obs  >= 0 & pa$prop_obs  <= 1 + 1e-10))
  expect_true(all(pa$prop_rand >= 0 & pa$prop_rand <= 1 + 1e-10))
})

test_that("fapa_pa is reproducible with same seed", {
  Xtilde <- make_ipsatized()
  pa1 <- fapa_pa(Xtilde, B = 100, seed = 7)
  pa2 <- fapa_pa(Xtilde, B = 100, seed = 7)
  expect_equal(pa1$n_retain, pa2$n_retain)
  expect_equal(pa1$thresh,   pa2$thresh)
})

## ---- Stage 2: fapa_procrustes -----------------------------------------------

test_that("fapa_procrustes returns correct matrix dimensions", {
  Xtilde <- make_ipsatized()
  pr     <- fapa_procrustes(Xtilde, K = 2, B = 30, seed = 1)
  expect_equal(dim(pr$angles_mat), c(30L, 2L))
})

test_that("fapa_procrustes angles are in [0, 90] degrees", {
  Xtilde <- make_ipsatized()
  pr     <- fapa_procrustes(Xtilde, K = 2, B = 50, seed = 1)
  expect_true(all(pr$angles_mat >= 0 - 1e-8))
  expect_true(all(pr$angles_mat <= 90 + 1e-8))
})

test_that("fapa_procrustes prop_stable is in [0, 1]", {
  Xtilde <- make_ipsatized(P = 200)
  pr     <- fapa_procrustes(Xtilde, K = 2, B = 50,
                             angle_thresh = 30, seed = 1)
  expect_gte(pr$prop_stable, 0)
  expect_lte(pr$prop_stable, 1)
})

## ---- Stage 3: fapa_tucker ---------------------------------------------------

test_that("fapa_tucker CC matrix dimensions are correct", {
  Xtilde <- make_ipsatized()
  tc     <- fapa_tucker(Xtilde, K = 2, B = 30, seed = 1)
  expect_equal(dim(tc$cc_mat), c(30L, 2L))
})

test_that("fapa_tucker CCs are in [-1, 1]", {
  Xtilde <- make_ipsatized()
  tc     <- fapa_tucker(Xtilde, K = 2, B = 50, seed = 1)
  expect_true(all(tc$cc_mat >= -1 - 1e-8, na.rm = TRUE))
  expect_true(all(tc$cc_mat <=  1 + 1e-8, na.rm = TRUE))
})

test_that("fapa_tucker mean CCs are high for stable data", {
  ## Key design: force V columns to be mean-zero so the 2-component signal
  ## survives row-centering (ipsatization) exactly.
  ## When colMeans(V) == 0, rowMeans(U %*% t(V)) == 0 for every person,
  ## so ipsatization leaves the signal unchanged and CCs are reliably high.
  set.seed(99)
  P <- 300;  I <- 10
  V <- matrix(rnorm(I * 2), I, 2)
  V <- scale(V, center = TRUE, scale = FALSE)   # force mean-zero columns
  U <- matrix(rnorm(P * 2), P, 2)
  signal <- U %*% t(V)
  noise  <- matrix(rnorm(P * I, sd = 0.15), P, I)
  Xtilde <- signal + noise
  Xtilde <- Xtilde - rowMeans(Xtilde)

  tc <- fapa_tucker(Xtilde, K = 2, B = 200, seed = 1)
  expect_true(all(tc$cc_mean > 0.90))
})
