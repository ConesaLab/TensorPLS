testthat::skip_if_not_installed("MASS")

library(testthat)


make_groups <- function(n) {
  set.seed(1)
  factor(sample(c("A","B"), n, replace = TRUE))
}

sim_X <- function(n, p, k, groups, effectX = 0.6, seed = 2) {
  set.seed(seed)
  X <- array(rnorm(n * p * k), dim = c(n, p, k))
  # small class effect on first 2 vars across first 2 slices
  v_eff <- seq_len(min(2, p))
  t_eff <- seq_len(min(2, k))
  for (t in t_eff) for (v in v_eff) X[groups=="B", v, t] <- X[groups=="B", v, t] + effectX
  dimnames(X) <- list(
    sample = paste0("S", seq_len(n)),
    variable = paste0("X", seq_len(p)),
    time = paste0("t", seq_len(k))
  )
  X
}

# Classical mode uses a NUMERIC matrix (not a factor/data.frame)
sim_Y_classic_num <- function(n, groups, q = 1, effectY = 0.8, seed = 3) {
  set.seed(seed)
  Y <- matrix(rnorm(n * q), n, q, dimnames = list(paste0("S", seq_len(n)), paste0("Y", seq_len(q))))
  Y[groups=="B", 1] <- Y[groups=="B", 1] + effectY
  Y
}

# Regression class-2: 3D numeric response with q channels and k slices
sim_Y_class2 <- function(n, q, k, groups, effectY = 0.5, seed = 4) {
  set.seed(seed)
  Y <- array(rnorm(n * q * k), dim = c(n, q, k))
  # inject class effect in first response on first slice
  Y[groups=="B", 1, 1] <- Y[groups=="B", 1, 1] + effectY
  dimnames(Y) <- list(
    sample = paste0("S", seq_len(n)),
    response = paste0("Y", seq_len(q)),
    time = paste0("t", seq_len(k))
  )
  Y
}

comp_names <- function(ncomp) paste0("Comp", seq_len(ncomp))

## ==========================================================
## 1) Classical mode: shapes + labels
## ==========================================================

test_that("Classical mode: factors shapes and labels are consistent", {
  n <- 30; p <- 6; k <- 3; q <- 1; ncomp <- 3
  groups <- make_groups(n)
  X <- sim_X(n, p, k, groups)
  Y <- sim_Y_classic_num(n, groups, q = q)
  
  res <- compute_npls_factors(X, Y, ncomp = ncomp)
  
  # FactorsX (Tt: n x ncomp; WsupraJ: p x ncomp; WsupraK: k x ncomp)
  expect_equal(dim(res$FactorsX$Mode1), c(n, ncomp))      # Tt
  expect_equal(dim(res$FactorsX$Mode2), c(p, ncomp))      # WsupraJ
  expect_equal(dim(res$FactorsX$Mode3), c(k, ncomp))      # WsupraK
  
  # FactorsY is NULL in classical mode
  expect_null(res$FactorsY)
  
  # Labels
  expect_identical(rownames(res$FactorsX$Mode1), dimnames(X)[[1]])
  expect_identical(rownames(res$FactorsX$Mode2), dimnames(X)[[2]])
  expect_identical(rownames(res$FactorsX$Mode3), dimnames(X)[[3]])
  expect_identical(colnames(res$FactorsX$Mode1), comp_names(ncomp))
  expect_identical(colnames(res$FactorsX$Mode2), comp_names(ncomp))
  expect_identical(colnames(res$FactorsX$Mode3), comp_names(ncomp))
  
  # B, G, Gu
  expect_equal(dim(res$B), c(ncomp, ncomp))
  expect_identical(rownames(res$B), comp_names(ncomp))
  expect_identical(colnames(res$B), comp_names(ncomp))
  expect_equal(dim(res$G), c(ncomp, ncomp, ncomp))
  expect_length(res$Gu, ncomp)
  expect_true(all(vapply(res$Gu, is.matrix, logical(1))))
})

## ==========================================================
## 2) Regression class-2: shapes + labels
## ==========================================================

test_that("Regression class-2: factors shapes and labels are consistent", {
  n <- 28; p <- 5; k <- 3; q <- 2; ncomp <- 3
  groups <- make_groups(n)
  X  <- sim_X(n, p, k, groups)
  Y3 <- sim_Y_class2(n, q, k, groups)
  
  res <- compute_npls_factors(X, Y3, ncomp = ncomp)
  
  # FactorsX
  expect_equal(dim(res$FactorsX$Mode1), c(n, ncomp))  # Tt
  expect_equal(dim(res$FactorsX$Mode2), c(p, ncomp))  # WsupraJ
  expect_equal(dim(res$FactorsX$Mode3), c(k, ncomp))  # WsupraK
  
  # FactorsY must exist in class-2
  expect_type(res$FactorsY, "list")
  expect_equal(dim(res$FactorsY$Mode1), c(n, ncomp))  # U
  expect_equal(dim(res$FactorsY$Mode2), c(q, ncomp))  # QsupraJ
  expect_equal(dim(res$FactorsY$Mode3), c(k, ncomp))  # QsupraK
  
  # Labels X
  expect_identical(rownames(res$FactorsX$Mode1), dimnames(X)[[1]])
  expect_identical(rownames(res$FactorsX$Mode2), dimnames(X)[[2]])
  expect_identical(rownames(res$FactorsX$Mode3), dimnames(X)[[3]])
  expect_identical(colnames(res$FactorsX$Mode1), comp_names(ncomp))
  expect_identical(colnames(res$FactorsX$Mode2), comp_names(ncomp))
  expect_identical(colnames(res$FactorsX$Mode3), comp_names(ncomp))
  
  # Labels Y
  expect_identical(rownames(res$FactorsY$Mode1), dimnames(X)[[1]])   # U rows = samples
  expect_identical(rownames(res$FactorsY$Mode2), dimnames(Y3)[[2]])  # responses
  expect_identical(rownames(res$FactorsY$Mode3), dimnames(Y3)[[3]])  # time slices
  expect_identical(colnames(res$FactorsY$Mode1), comp_names(ncomp))
  expect_identical(colnames(res$FactorsY$Mode2), comp_names(ncomp))
  expect_identical(colnames(res$FactorsY$Mode3), comp_names(ncomp))
  
  # B, G, Gu
  expect_equal(dim(res$B), c(ncomp, ncomp))
  expect_equal(dim(res$G), c(ncomp, ncomp, ncomp))
  expect_length(res$Gu, ncomp)
})

## ==========================================================
## 3) Validation errors (bounds, feasibility, structure, NA)
## ==========================================================

test_that("Input validation: bounds, feasibility, structure, NA", {
  n <- 12; p <- 4; k <- 2; q <- 2
  groups <- make_groups(n)
  X  <- sim_X(n, p, k, groups)
  Y2 <- sim_Y_classic_num(n, groups, q = 1)
  Y3 <- sim_Y_class2(n, q, k, groups)
  
  # ncomp out of allowed [2..10]
  expect_error(
    compute_npls_factors(X, Y3, ncomp = 1),
    "between 2 and 10", fixed = TRUE
  )
  expect_error(
    compute_npls_factors(X, Y3, ncomp = 11),
    "between 2 and 10", fixed = TRUE
  )
  
  # Statistical feasibility: min(n-1, p) = min(11,4) = 4 -> choose 5
  expect_error(
    compute_npls_factors(X, Y2, ncomp = 5),
    "exceeds the statistical limit", fixed = TRUE
  )
  
  # X must be 3D
  expect_error(
    compute_npls_factors(X = X[, , 1], Y = Y2, ncomp = 2),
    "X must be a 3D array", ignore.case = TRUE
  )
  
  # n mismatch
  expect_error(
    compute_npls_factors(X[1:10, , , drop = FALSE], Y2, ncomp = 2),
    "same number of observations", ignore.case = TRUE
  )
  
  # NA guard
  Xna <- X
  Xna[1, 1, 1] <- NA
  expect_error(
    compute_npls_factors(Xna, Y3, ncomp = 2),
    "must not contain NA", ignore.case = TRUE
  )
})

## ==========================================================
## 4) Minimal smoke test on Gu/G contents (finite) 
## ==========================================================

test_that("Gu and G contain finite values and expected lengths", {
  n <- 20; p <- 5; k <- 3; q <- 2; ncomp <- 2
  groups <- make_groups(n)
  X  <- sim_X(n, p, k, groups)
  Y3 <- sim_Y_class2(n, q, k, groups)
  
  res <- compute_npls_factors(X, Y3, ncomp = ncomp)
  
  expect_length(res$Gu, ncomp)
  expect_true(all(vapply(res$Gu, function(M) all(is.finite(M)), logical(1))))
  expect_true(all(is.finite(res$G)))
  expect_identical(dimnames(res$G),
                   rep(list(comp_names(ncomp)), 3))
})
