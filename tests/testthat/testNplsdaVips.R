
context("Testing nplsda-core")

skip_if_not_installed("mixOmics")

set.seed(42)

# --------- Synthetic data generators ---------

gen_tensor_classic <- function(n = 40, p = 6, k = 3, qY = 2, noise_sd = 0.3,
                               with_dimnames = FALSE) {
  # X: n x p x k
  X <- array(rnorm(n * p * k), dim = c(n, p, k))
  # Build Y as linear combination across slices + noise
  Y <- matrix(0, n, qY)
  for (tt in seq_len(k)) {
    B <- matrix(runif(p * qY, -1, 1), p, qY)
    Y <- Y + X[, , tt] %*% B
  }
  Y <- Y + matrix(rnorm(n * qY, sd = noise_sd), n, qY)
  
  if (with_dimnames) {
    dimnames(X) <- list(NULL,
                        paste0("Gene", seq_len(p)),
                        paste0("T", seq_len(k)))
  }
  list(X = X, Y = Y)
}

gen_tensor_class2 <- function(n = 32, p = 5, k = 3, q = 2, noise_sd = 0.25,
                              with_dimnames = FALSE) {
  # X: n x p x k
  X <- array(rnorm(n * p * k), dim = c(n, p, k))
  # Y(3D): n x q x k
  Y <- array(0, dim = c(n, q, k))
  for (tt in seq_len(k)) {
    B <- matrix(runif(p * q, -1, 1), p, q)
    Y[, , tt] <- X[, , tt] %*% B + matrix(rnorm(n * q, sd = noise_sd), n, q)
  }
  
  if (with_dimnames) {
    dimnames(X) <- list(NULL,
                        paste0("Var", seq_len(p)),
                        paste0("Time", seq_len(k)))
  }
  list(X = X, Y = Y)
}


# --------- Tests ---------

test_that("Input validation errors are informative", {
  # What this test checks:
  # - X not 3D -> error
  # - row mismatch between X and Y -> error
  # - ncomp outside [2, 10] -> error
  # - ncomp > min(n-1, p) -> error with statistical limit message
  # - NA presence -> error
  # - outcome.Y provided but Y not 3D -> error
  # - Y is 3D with q>1 and outcome.Y missing -> error from main function
  
  dat <- gen_tensor_classic()
  X3 <- dat$X; Ymat <- dat$Y
  
  # X not 3D
  expect_error(nplsda_vips(X = matrix(1, 2, 2), Y = Ymat),
               "X should be a 3D array")
  
  # Row mismatch
  expect_error(nplsda_vips(X = X3[1:30, , , drop = FALSE], Y = Ymat),
               "X and Y must have the same number of observations")
  
  # ncomp outside [2, 10]
  expect_error(nplsda_vips(X3, Ymat, ncomp = 1),
               "ncomp must be between 2 and 10")
  
  # ncomp > min(n-1, p)
  n <- dim(X3)[1]; p <- dim(X3)[2]
  too_big <- min(n - 1, p) + 1
  expect_error(nplsda_vips(X3, Ymat, ncomp = too_big),
               "exceeds the statistical limit")
  
 
  
  # NA presence
  X_bad <- X3; X_bad[1, 1, 1] <- NA_real_
  expect_error(nplsda_vips(X_bad, Ymat),
               "cannot contain NA values")
  
  # outcome.Y provided but Y not 3D
  expect_error(nplsda_vips(X3, Ymat, outcome.Y = factor(rep(1, nrow(Ymat)))),
               "outcome.Y was provided but Y is not 3D")
  
  # Y is 3D with q>1 and outcome.Y missing
  dat2 <- gen_tensor_class2()
  expect_error(nplsda_vips(dat2$X, dat2$Y),
               "outcome.Y required for regression class 2")
})

test_that("Classic mode returns consistent structure, shapes, and finite values", {
  # What this test checks:
  # - S3 class 'nplsda_vips'
  # - Core fields exist with expected dimensions
  # - VIP matrices non-negative & finite
  # - Q2/Q2cum shapes consistent with Y
  # - y.pred/residuals finite and correctly shaped
  # - regression_class2 == FALSE
  
  dat <- gen_tensor_classic(n = 36, p = 5, k = 3, qY = 2, noise_sd = 0.2)
  X <- dat$X; Y <- dat$Y
  ncomp <- 3
  
  res <- nplsda_vips(X, Y, ncomp = ncomp, slice_vip = FALSE)
  
  expect_s3_class(res, "nplsda_vips")
  
  # VIP2D: (p*k) x ncomp
  pk <- dim(X)[2] * dim(X)[3]
  expect_equal(dim(res$VIP2D), c(pk, ncomp))
  expect_equal(colnames(res$VIP2D), paste0("Comp", seq_len(ncomp)))
  
  # VIP3Dmodel1: p x ncomp ; VIP3Dmodel2: p x k
  expect_equal(dim(res$VIP3Dmodel1), c(dim(X)[2], ncomp))
  expect_equal(dim(res$VIP3Dmodel2), c(dim(X)[2], dim(X)[3]))
  
  # VIP non-negative & finite
  expect_true(all(is.finite(res$VIP2D)))
  expect_true(all(res$VIP2D >= 0))
  expect_true(all(is.finite(res$VIP3Dmodel1)))
  expect_true(all(res$VIP3Dmodel1 >= 0))
  expect_true(all(is.finite(res$VIP3Dmodel2)))
  expect_true(all(res$VIP3Dmodel2 >= 0))
  
  # Q2 / Q2cum: ncomp x q
  expect_equal(dim(res$Q2),    c(ncomp, ncol(Y)))
  expect_equal(dim(res$Q2cum), c(ncomp, ncol(Y)))
  expect_true(all(is.finite(res$Q2)))
  expect_true(all(is.finite(res$Q2cum)))
  expect_true(all(res$Q2 <= 1 + 1e-8))
  
  # Mean Q2 (3D) metrics: ncomp x 1
  expect_equal(dim(res$NPLSDAQ2mean3D), c(ncomp, 1))
  expect_equal(dim(res$NPLSDAQ2cummean3D), c(ncomp, 1))
  
  # Explained variance shapes and ranges
  expect_equal(dim(res$Explvar3D), c(ncomp, 2))
  expect_equal(colnames(res$Explvar3D), c("R2Xmean3D", "R2Ymean3D"))
  expect_true(all(res$Explvar3D >= -1e-8))
  expect_true(all(res$Explvar3D <= 1 + 1e-8))
  
  expect_equal(dim(res$NPLSDAexplVar), c(ncomp, 1))
  expect_equal(colnames(res$NPLSDAexplVar), "R2.Y")
  
  expect_equal(dim(res$explvar), c(ncomp, 4))
  expect_equal(colnames(res$explvar), c("R2X", "R2Xcum", "R2Y", "R2Ycum"),
               ignore_attr = TRUE) # names could differ depending on code edits
  
  # Predictions & residuals
  expect_equal(dim(res$y.pred), c(dim(X)[1], ncol(Y)))
  expect_equal(dim(res$residuals), c(dim(X)[1], ncol(Y)))
  expect_true(all(is.finite(res$y.pred)))
  expect_true(all(is.finite(res$residuals)))
  
  expect_false(res$regression_class2)
})

test_that("VIP naming uses X dimnames (var_time) when provided", {
  # What this test checks:
  # - If X has names for mode-2 and mode-3, VIP2D rownames are "Var_Time".
  # - VIP3Dmodel2 colnames match the 3rd-mode names.
  dat <- gen_tensor_classic(with_dimnames = TRUE)
  X <- dat$X; Y <- dat$Y
  ncomp <- 3
  
  res <- nplsda_vips(X, Y, ncomp = ncomp, slice_vip = FALSE)
  
  var_names  <- dimnames(X)[[2]]
  time_names <- dimnames(X)[[3]]
  
  expected_rownames <- as.vector(outer(var_names, time_names, paste, sep = "_"))
  expect_identical(rownames(res$VIP2D), expected_rownames)
  expect_identical(colnames(res$VIP3Dmodel2), time_names)
  expect_equal(colnames(res$VIP3Dmodel1), paste0("t", seq_len(ncomp), "mean3D"))
})

test_that("slice_vip = TRUE computes per-slice VIP aggregates with correct shapes", {
  # What this test checks:
  # - VIP3Dmodel1 still p x ncomp, VIP3Dmodel2 still p x k
  # - Values remain finite and non-negative
  dat <- gen_tensor_classic(n = 36, p = 5, k = 4, qY = 2, noise_sd = 0.25)
  X <- dat$X; Y <- dat$Y
  ncomp <- 3
  
  res <- nplsda_vips(X, Y, ncomp = ncomp,  slice_vip = TRUE)
  
  expect_equal(dim(res$VIP3Dmodel1), c(dim(X)[2], ncomp))
  expect_equal(dim(res$VIP3Dmodel2), c(dim(X)[2], dim(X)[3]))
  
  expect_true(all(is.finite(res$VIP3Dmodel1)))
  expect_true(all(res$VIP3Dmodel1 >= 0))
  expect_true(all(is.finite(res$VIP3Dmodel2)))
  expect_true(all(res$VIP3Dmodel2 >= 0))
})

##Note: if i test this with p < ncomp i have errors of .validate.
test_that("When k < ncomp, VIP3Dmodel1 still returns p x ncomp", {
  # What this test checks:
  # - The function uses slice_idx = min(ncomp, k) but returns a p x ncomp matrix.
  dat <- gen_tensor_classic(n = 30, p = 6, k = 2, qY = 1, noise_sd = 0.25)
  X <- dat$X; Y <- dat$Y
  ncomp <- 5
  
  res <- nplsda_vips(X, Y, ncomp = ncomp, slice_vip = FALSE)
  expect_equal(dim(res$VIP3Dmodel1), c(dim(X)[2], ncomp))
  expect_equal(colnames(res$VIP3Dmodel1), paste0("t", seq_len(ncomp), "mean3D"))
})

test_that("Regression class-2 mode works and flags regression_class2 = TRUE", {
  # What this test checks:
  # - With Y 3D and q>1 plus outcome.Y provided, function runs
  # - Shapes adapt: Q2 and y.pred use flattened Y (q*k columns)
  dat <- gen_tensor_class2(n = 28, p = 4, k = 3, q = 2, noise_sd = 0.25)
  X <- dat$X; Y3 <- dat$Y
  ncomp <- 3
  
  # outcome for class-2 (arbitrary labels)
  outcome <- factor(sample(letters[1:2], size = dim(X)[1], replace = TRUE))
  
  res <- nplsda_vips(X, Y3, ncomp = ncomp, slice_vip = FALSE,
                     outcome.Y = outcome)
  
  expect_true(res$regression_class2)
  
  q <- dim(Y3)[2]; k <- dim(Y3)[3]
  expect_equal(dim(res$Q2), c(ncomp, q * k))
  expect_equal(dim(res$Q2cum), c(ncomp, q * k))
  expect_equal(dim(res$y.pred), c(dim(X)[1], q * k))
})

test_that("Deterministic output with fixed RNG seed (CV fold assignment)", {
  # What this test checks:
  # - With a fixed seed, CV-based quantities (e.g., Q2) are reproducible.
  dat <- gen_tensor_classic(n = 34, p = 5, k = 3, qY = 2, noise_sd = 0.2)
  X <- dat$X; Y <- dat$Y
  ncomp <- 3
  
  set.seed(777)
  a <- nplsda_vips(X, Y, ncomp = ncomp, slice_vip = FALSE)
  
  set.seed(777)
  b <- nplsda_vips(X, Y, ncomp = ncomp, slice_vip = FALSE)
  
  expect_equal(a$Q2,    b$Q2,    tolerance = 1e-12)
  expect_equal(a$Q2cum, b$Q2cum, tolerance = 1e-12)
  expect_equal(a$NPLSDAQ2mean3D, b$NPLSDAQ2mean3D, tolerance = 1e-12)
  expect_equal(a$NPLSDAQ2cummean3D, b$NPLSDAQ2cummean3D, tolerance = 1e-12)
})

