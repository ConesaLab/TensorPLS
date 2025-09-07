
context("validates the core behavior, edge cases, and invariants of `pls_reg")

set.seed(123)

gen_data <- function(n = 60, p = 8, q = 2, noise_sd = 0.3) {
  X <- matrix(rnorm(n * p), n, p)
  B <- matrix(runif(p * q, -1, 1), p, q)
  Y_mean <- matrix(rnorm(q), 1, q)
  Y <- X %*% B + matrix(rnorm(n * q, sd = noise_sd), n, q) + 
    matrix(rep(Y_mean, each = n), n, q)
  list(X = X, Y = Y, B = B)
}

test_that("Input validation errors are informative", {
  # checks:
  # 1) Non-matrix inputs raise an error.
  # 2) NA in inputs raises an error.
  # 3) Row mismatch between X and Y raises an error.
  # 4) ncomp < 2 raises an error.
  
  dat <- gen_data()
  X <- dat$X; Y <- dat$Y
  
  expect_error(pls_reg(as.data.frame(X), Y), "X should be a numeric matrix")
  X_bad <- X; X_bad[1,1] <- NA_real_
  expect_error(pls_reg(X_bad, Y), "X contain NA")
  expect_error(pls_reg(X, Y[-1, , drop = FALSE]), "X and Y: #differents rows")
  expect_error(pls_reg(X, Y, ncomp = 1), "ncomp should be ≥2")
})

test_that("Basic shape, class, and component matrices are consistent (cv = TRUE)", {
  # checks:
  # - Matrices have expected dimensions.
  # - No NA/Inf in core outputs when data carry signal.
  # - Q2 and Q2cum are present and have expected shape when cv = TRUE.
  
  dat <- gen_data(n = 50, p = 6, q = 2, noise_sd = 0.2)
  X <- dat$X; Y <- dat$Y
  obj <- pls_reg(X, Y, ncomp = 4, cv = TRUE)
  

  nc <- 4
  expect_equal(dim(obj$x.scores), c(nrow(X), nc))
  expect_equal(dim(obj$x.loads),  c(ncol(X), nc))
  expect_equal(dim(obj$y.scores), c(nrow(X), nc))
  expect_equal(dim(obj$y.loads),  c(ncol(Y), nc))
  expect_equal(dim(obj$raw.wgs),  c(ncol(X), nc))
  expect_equal(dim(obj$mod.wgs),  c(ncol(X), nc))
  
  expect_equal(dim(obj$std.coefs), c(ncol(X), ncol(Y)))
  expect_equal(dim(obj$reg.coefs), c(ncol(X) + 1, ncol(Y)))  # + intercept row
  
  expect_equal(dim(obj$y.pred), c(nrow(X), ncol(Y)))
  expect_equal(dim(obj$resid),  c(nrow(X), ncol(Y)))
  
  expect_equal(dim(obj$VIP), c(ncol(X), nc))
  
  # CV-specific components
  expect_equal(dim(obj$Q2),    c(nc, ncol(Y)))
  expect_equal(dim(obj$Q2cum), c(nc, ncol(Y)))
  
  # Sanity: finite values in key outputs
  expect_true(all(is.finite(obj$y.pred)))
  expect_true(all(is.finite(obj$VIP)))
  expect_true(all(is.finite(obj$Q2)))
  expect_true(all(is.finite(obj$Q2cum)))
})

test_that("Residual + prediction equals original Y (within tolerance)", {
  #  checks:
  # y.pred + resid must reconstruct Y (definition of residuals).
  
  dat <- gen_data()
  X <- dat$X; Y <- dat$Y
  obj <- pls_reg(X, Y, ncomp = 3, cv = TRUE)
  recon <- obj$y.pred + obj$resid
  expect_equal(recon, Y, tolerance = 1e-7)
})

test_that("Q2cum matches the multiplicative definition from Q2", {
  # W checks:
  # Q2cum[h, ] == 1 - prod_{j=1..h} (1 - Q2[j, ])
  
  dat <- gen_data()
  X <- dat$X; Y <- dat$Y
  obj <- pls_reg(X, Y, ncomp = 3, cv = TRUE)
  
  Q2 <- obj$Q2
  Q2cum_expected <- apply(Q2, 2, function(col) {
    sapply(seq_along(col), function(h) 1 - prod(1 - col[1:h]))
  })
  # Ensure shape (nc x q)
  if (is.null(dim(Q2cum_expected))) {
    Q2cum_expected <- matrix(Q2cum_expected, nrow = nrow(Q2), ncol = ncol(Q2))
  }
  expect_equal(obj$Q2cum, Q2cum_expected, tolerance = 1e-7)
})

test_that("VIP matrix is non-negative and finite under a clear signal", {
  #  checks:
  # VIPs should be finite and non-negative ().
  
  dat <- gen_data(noise_sd = 0.15)
  X <- dat$X; Y <- dat$Y
  obj <- pls_reg(X, Y, ncomp = 5, cv = TRUE)
  
  expect_true(all(is.finite(obj$VIP)))
  expect_true(all(obj$VIP >= 0))
})

test_that("Component count is truncated to min(n-1, p) when ncomp is too large", {
  # checks:
  # If requested ncomp exceeds n-1 or p, the function uses nc = min(ncomp, n-1, p).
  
  n <- 30; p <- 5; q <- 2
  dat <- gen_data(n = n, p = p, q = q)
  X <- dat$X; Y <- dat$Y
  
  obj <- pls_reg(X, Y, ncomp = 20, cv = TRUE)
  # Expected nc:
  nc <- min(20, n - 1, p)
  
  expect_equal(ncol(obj$x.scores), nc)
  expect_equal(ncol(obj$x.loads),  nc)
  expect_equal(ncol(obj$y.scores), nc)
  expect_equal(ncol(obj$y.loads),  nc)
  expect_equal(ncol(obj$raw.wgs),  nc)
  expect_equal(ncol(obj$mod.wgs),  nc)
  expect_equal(ncol(obj$VIP),      nc)
  expect_equal(nrow(obj$Q2),       nc)
  expect_equal(nrow(obj$Q2cum),    nc)
})

test_that("cv = FALSE disables Q2/Q2cum and still returns consistent components", {
  # checks:
  # With cv = FALSE, Q2 and Q2cum are NULL; other outputs are still valid.
  
  dat <- gen_data()
  X <- dat$X; Y <- dat$Y
  obj <- pls_reg(X, Y, ncomp = 3, cv = FALSE)
  
  expect_null(obj$Q2)
  expect_null(obj$Q2cum)
  
  expect_equal(dim(obj$x.scores), c(nrow(X), 3))
  expect_equal(dim(obj$VIP), c(ncol(X), 3))
  expect_true(all(is.finite(obj$y.pred)))
})

test_that("Non-logical cv input is treated as FALSE", {
  #  checks:
  # Passing a non-logical cv (e.g., integer 1) coerces to FALSE per implementation,
  # thus Q2 and Q2cum should be NULL.
  
  dat <- gen_data()
  X <- dat$X; Y <- dat$Y
  obj <- pls_reg(X, Y, ncomp = 3, cv = 1L)  # not a logical
  expect_null(obj$Q2)
  expect_null(obj$Q2cum)
})

test_that("Results are deterministic given a fixed seed (k-fold sampling)", {
  #checks:
  # The internal fold assignment uses random sampling. With a fixed seed,
  # repeated runs should produce identical CV-based metrics.
  
  dat <- gen_data()
  X <- dat$X; Y <- dat$Y
  
  set.seed(999)
  a <- pls_reg(X, Y, ncomp = 4, cv = TRUE)
  
  set.seed(999)
  b <- pls_reg(X, Y, ncomp = 4, cv = TRUE)
  
  expect_equal(a$Q2,    b$Q2,    tolerance = 1e-12)
  expect_equal(a$Q2cum, b$Q2cum, tolerance = 1e-12)
})

