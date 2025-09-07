library(testthat)

# ---- tiny helpers
make_lowrank <- function(I=6, J=5, K=4, P=3, Q=2, R=2, noise_sd=0) {
  set.seed(1)
  qrQ <- function(M) qr.Q(qr(M))  # orthonormal columns
  A <- qrQ(matrix(rnorm(I*P), I, P))
  B <- qrQ(matrix(rnorm(J*Q), J, Q))
  C <- qrQ(matrix(rnorm(K*R), K, R))
  G <- array(rnorm(P*Q*R), c(P,Q,R))
  
  # reconstruct X = G x1 A x2 B x3 C
  X <- array(0, c(I,J,K))
  for (i in 1:I) for (j in 1:J) for (k in 1:K) {
    s <- 0
    for (p in 1:P) for (q in 1:Q) for (r in 1:R)
      s <- s + A[i,p]*B[j,q]*C[k,r]*G[p,q,r]
    X[i,j,k] <- s
  }
  if (noise_sd > 0) X <- X + array(rnorm(I*J*K, sd=noise_sd), c(I,J,K))
  list(X=X, A=A, B=B, C=C, G=G)
}

# ---- Imputemethod

test_that("Imputemethod: input checks and no-NA short-circuit", {
  X <- array(rnorm(6*5*4), c(6,5,4))
  expect_no_condition(Imputemethod(X, verbose=FALSE))
  # no-NA -> identical object
  out <- Imputemethod(X, verbose=FALSE)
  expect_identical(out, X)
})

test_that("Imputemethod: imputes only missing entries; observed stay unchanged; seed reproducible", {
  withr::local_seed(42)
  toy <- make_lowrank()$X
  idx <- sample(length(toy), size = 0.1*length(toy))
  Xna <- toy; Xna[idx] <- NA
  
  out1 <- Imputemethod(Xna, fac=c(3,2,2), conver=1e-6, max.iter=200, seed=123, verbose=FALSE)
  out2 <- Imputemethod(Xna, fac=c(3,2,2), conver=1e-6, max.iter=200, seed=123, verbose=FALSE)
  
  # no missing left
  expect_false(anyNA(out1))
  # observed entries preserved
  expect_equal(out1[-idx], toy[-idx], tolerance = 0)  # exact equality
  # reproducible with same seed
  expect_equal(out1, out2, tolerance = 0)
})

test_that("Imputemethod: fac clipped and sd=0 edge-case handled", {
  Xconst <- array(5, c(4,3,2)); Xconst[1,1,1] <- NA
  expect_no_condition(Imputemethod(Xconst, fac=c(9,9,9), seed=1, verbose=FALSE))
  out <- Imputemethod(Xconst, fac=c(9,9,9), seed=1, verbose=FALSE)
  expect_false(anyNA(out))
})

# ---- tucker3mod2-

test_that("tucker3mod2: errors if COMP exceeds dims; basic input checks", {
  X <- array(rnorm(6*5*4), c(6,5,4))
  expect_error(tucker3mod2(X, COMP=c(10,2,2)), "exceed")
  expect_error(tucker3mod2(matrix(1, 2, 2), COMP=c(1,1,1)))  # not a 3-way array
})

test_that("tucker3mod2: internal imputation path triggers with NA (message)", {
  toy <- make_lowrank()$X
  toy[1,1,1] <- NA
  expect_message(
    tucker3mod2(toy, COMP=c(3,2,2), seed=1, verbose=TRUE),
    "Missing values detected"
  )
})

test_that("tucker3mod2: near-perfect reconstruction on noiseless low-rank", {
  testthat::skip_on_cran()
  
  toy <- make_lowrank(I=6,J=5,K=4,P=3,Q=2,R=2, noise_sd=0)
  fit <- tucker3mod2(toy$X, COMP=c(3,2,2), conver=1e-8, max.iter=2000, seed=7, verbose=FALSE)
  
  # identities and high explained variance
  expect_lt(fit$SSE, 1e-8 * fit$SST)
  expect_gt(fit$expl.var, 0.999999)
  expect_equal(fit$SSF, sum(fit$Xhat^2), tolerance = 1e-8)
  expect_equal(fit$SST, fit$SSF + fit$SSE, tolerance = 1e-8)
  
  # factors roughly orthonormal
  expect_equal(crossprod(fit$FactorsX$Mode1), diag(3), tolerance = 1e-6)
  expect_equal(crossprod(fit$FactorsX$Mode2), diag(2), tolerance = 1e-6)
  expect_equal(crossprod(fit$FactorsX$Mode3), diag(2), tolerance = 1e-6)
})

test_that("tucker3mod2: reproducible with seed", {
  toy <- make_lowrank()$X
  f1 <- tucker3mod2(toy, COMP=c(3,2,2), seed=123, max.iter=300, verbose=FALSE)
  f2 <- tucker3mod2(toy, COMP=c(3,2,2), seed=123, max.iter=300, verbose=FALSE)
  expect_equal(f1$expl.var, f2$expl.var, tolerance = 0)  # identical with same seed
})

