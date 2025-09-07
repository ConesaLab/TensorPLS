
# ---------------------------------------------------------------

testthat::skip_if_not_installed("mixOmics")
testthat::skip_if_not_installed("abind")

library(testthat)


make_groups <- function(n) {
  factor(sample(c("A","B"), n, replace = TRUE))
}

sim_X <- function(n, p, k, groups, effectX = 0.6, seed = 1) {
  set.seed(seed)
  X <- array(rnorm(n * p * k, sd = 1), dim = c(n, p, k))
  # Inject a small class effect on first 2 variables and first 2 slices
  idx_v <- seq_len(min(2, p))
  idx_t <- seq_len(min(2, k))
  for (t in idx_t) {
    for (v in idx_v) {
      X[groups == "B", v, t] <- X[groups == "B", v, t] + effectX
    }
  }
  # Attach dimnames
  dimnames(X) <- list(
    sample = paste0("S", seq_len(n)),
    variable = paste0("X", seq_len(p)),
    time = paste0("t", seq_len(k))
  )
  X
}

# Y for classical mode: data.frame so plsda sees a factor in first column
sim_Y_classic <- function(groups) {
  data.frame(group = groups)
}

# Y for regression class-2 (3D): (n x q x k) with dimnames
sim_Y_class2 <- function(n, q, k, groups, effectY = 0.5, seed = 2) {
  set.seed(seed)
  Y <- array(rnorm(n * q * k, sd = 1), dim = c(n, q, k))
  # Inject class effect in first response channel on first slice
  q_eff <- min(1, q)
  k_eff <- min(1, k)
  Y[groups == "B", seq_len(q_eff), seq_len(k_eff)] <-
    Y[groups == "B", seq_len(q_eff), seq_len(k_eff)] + effectY
  dimnames(Y) <- list(
    sample = paste0("S", seq_len(n)),
    response = paste0("Y", seq_len(q)),
    time = paste0("t", seq_len(k))
  )
  Y
}

#: expect a matrix with at least ncomp columns
# matrix has at least ncomp columns (use for GLOBAL checks)
.expect_matrix_mincomp <- function(M, ncomp) {
  expect_true(is.matrix(M))
  expect_gte(ncol(M), as.integer(ncomp))
}

# matrix has exactly the slice-wise effective ncomp (use for CLASS-2 SLICE checks)
.expect_slice_block_ncomp <- function(M, n, p, q, ncomp) {
  expect_true(is.matrix(M))
  eff <- min(as.integer(ncomp), n - 1L, p, q)
  expect_equal(ncol(M), eff)
}

#  tests: classical mode

test_that("Classical mode: global, slice-wise, consensus shapes & labels", {
  n <- 30; p <- 6; k <- 3; ncomp <- 3
  groups <- make_groups(n)
  X <- sim_X(n, p, k, groups, effectX = 0.7)
  Y2 <- sim_Y_classic(groups)  # data.frame with factor in first column
  
  res <- compute_npls_loadings(X, Y2, ncomp = ncomp)
  
  # basic structure
  expect_s3_class(res, "npls_loadings_unified")
  expect_true(is.list(res$NPLSDAloadings))
  expect_true(all(c("X", "Y") %in% names(res$NPLSDAloadings)))
  
  # global loadings: rows = p*k for X; expect >= ncomp columns
  glX <- res$NPLSDAloadings$X
  .expect_matrix_mincomp(glX, ncomp)
  expect_equal(nrow(glX), p * k)
  
  # dimnames alignment: "var__t" labels should exist
  var_names  <- dimnames(X)[[2]]
  time_names <- dimnames(X)[[3]]
  expect_true(any(rownames(glX) %in% paste0(var_names[1], "__t", time_names[1])))
  
  # slice-wise (first slice)
  sw1 <- res$NPLSDAloadingsperMode3[["item:1"]][[1]]
  expect_true(is.list(sw1))
  expect_true(all(c("X","Y") %in% names(sw1)))
  .expect_matrix_mincomp(sw1$X, ncomp)
  expect_equal(nrow(sw1$X), p)  # per-slice X loadings have p rows
  
  # consensus: rows = p; cols <= ncomp
  cons <- res$NPLSDAConsensusloadings
  expect_true(is.list(cons))
  expect_true(!is.null(cons$Block.XConsensus))
  expect_equal(nrow(cons$Block.XConsensus), p)
  expect_lte(ncol(cons$Block.XConsensus), ncomp)
  
  # consensus X rownames should match variable dimnames
  expect_true(all(rownames(cons$Block.XConsensus) == var_names))
  
  # transformed loadings (if present)
  if (!is.null(res$transformedNPLSDAloadings)) {
    expect_true(all(c("names","x","y") %in% colnames(res$transformedNPLSDAloadings)))
  }
})

# --- tests: regression class-2 mode

test_that("Regression class-2 mode: block loadings, slice-wise, consensus, labels", {
  n <- 30; p <- 5; q <- 2; k <- 3; ncomp <- 3
  groups <- make_groups(n)
  X <- sim_X(n, p, k, groups, effectX = 0.6)
  Y3 <- sim_Y_class2(n, q, k, groups, effectY = 0.5)
  
  res <- compute_npls_loadings(X, Y3, ncomp = ncomp, outcome.Y = groups)
  
  # global loadings: block structure
  expect_true(is.list(res$NPLSDAloadings))
  nm <- names(res$NPLSDAloadings)
  expect_true(all(c("Block.X","Block.Y") %in% nm))
  
  glBX <- res$NPLSDAloadings[["Block.X"]]
  glBY <- res$NPLSDAloadings[["Block.Y"]]
  
  .expect_matrix_mincomp(glBX, ncomp)
  .expect_matrix_mincomp(glBY, ncomp)
  expect_equal(nrow(glBX), p * k)
  expect_equal(nrow(glBY), q * k)
  
  # dimnames alignment for unfolded labels
  var_names  <- dimnames(X)[[2]]
  time_names <- dimnames(X)[[3]]
  y_names    <- dimnames(Y3)[[2]]
  y_times    <- dimnames(Y3)[[3]]
  
  expect_true(any(rownames(glBX) %in% paste0(var_names[1], "__t", time_names[1])))
  expect_true(any(rownames(glBY) %in% paste0(y_names[1], "__t", y_times[1])))
  
  # slice-wise (first slice): per-block matrices with p and q rows
  sw1 <- res$NPLSDAloadingsperMode3[["item:1"]][[1]]
  expect_true(is.list(sw1))
  expect_true(all(c("Block.X","Block.Y") %in% names(sw1)))
  expect_equal(nrow(sw1[["Block.X"]]), p)
  expect_equal(nrow(sw1[["Block.Y"]]), q)
  .expect_slice_block_ncomp(sw1[["Block.X"]], n, p, q, ncomp)
  .expect_slice_block_ncomp(sw1[["Block.Y"]], n, p, q, ncomp)
  # consensus: X has p rows w/ variable names; Y has q rows (row names may be NULL)
  cons <- res$NPLSDAConsensusloadings
  expect_equal(nrow(cons$Block.XConsensus), p)
  expect_lte(ncol(cons$Block.XConsensus), ncomp)
  expect_true(all(rownames(cons$Block.XConsensus) == var_names))
  
  expect_equal(nrow(cons$Block.YConsensus), q)
  expect_lte(ncol(cons$Block.YConsensus), ncomp)
})

# --- tests: error handling -------------------------------------

test_that("Errors on invalid inputs and missing outcome", {
  n <- 12; p <- 4; k <- 2; q <- 2
  groups <- make_groups(n)
  X <- sim_X(n, p, k, groups)
  Y2 <- sim_Y_classic(groups)                 # classical Y
  Y3 <- sim_Y_class2(n, q, k, groups)         # <-- define Y3 (3D) for class-2
  
  # Classical mode: ncomp above statistical limit but within [2,10]
  # max_stat = min(n-1, p) = min(11, 4) = 4  -> choose 5
  max_stat <- min(n - 1L, p)
  ncomp_bad <- min(10L, max_stat + 1L)        # = 5 here
  expect_error(
    compute_npls_loadings(X, Y2, ncomp = ncomp_bad),
    "exceeds the statistical limit",
    fixed = TRUE
  )
  
  # Regression class-2: ncomp out of bounds (too small)
  expect_error(
    compute_npls_loadings(X, Y3, ncomp = 1, outcome.Y = groups),
    "between 2 and 10", fixed = TRUE
  )
  
  # Regression class-2: missing outcome.Y
  expect_error(
    compute_npls_loadings(X, Y3, ncomp = 2),
    "outcome.Y.*required", ignore.case = TRUE
  )
})

# --- tests: dimnames pass-through on slices ---------------------

test_that("Slice-wise loadings keep variable names aligned (classical & class-2)", {
  n <- 20; p <- 5; k <- 3; q <- 2; ncomp <- 2
  groups <- make_groups(n)
  X <- sim_X(n, p, k, groups)
  
  # classical
  Y2 <- sim_Y_classic(groups)
  resC <- compute_npls_loadings(X, Y2, ncomp = ncomp)
  swC1 <- resC$NPLSDAloadingsperMode3[["item:1"]][[1]]
  expect_true(all(rownames(swC1$X) == dimnames(X)[[2]]))
  
  # regression class-2
  Y3 <- sim_Y_class2(n, q, k, groups)
  resR <- compute_npls_loadings(X, Y3, ncomp = ncomp, outcome.Y = groups)
  sw <- resR$NPLSDAloadingsperMode3[["item:1"]][[1]]
  .expect_slice_block_ncomp(sw[["Block.X"]], n, p, q, ncomp)
  .expect_slice_block_ncomp(sw[["Block.Y"]], n, p, q, ncomp)
})

test_that("Slice-wise loadings keep variable names aligned (classical & class-2)", {
  n <- 20; p <- 5; k <- 3; q <- 2; ncomp <- 2
  groups <- make_groups(n)
  X <- sim_X(n, p, k, groups)
  
  # classical
  Y2 <- sim_Y_classic(groups)
  resC <- compute_npls_loadings(X, Y2, ncomp = ncomp)
  swC1 <- resC$NPLSDAloadingsperMode3[["item:1"]][[1]]
  expect_true(all(rownames(swC1$X) == dimnames(X)[[2]]))
  
  # regression class-2
  Y3 <- sim_Y_class2(n, q, k, groups)
  resR <- compute_npls_loadings(X, Y3, ncomp = ncomp, outcome.Y = groups)
  sw1 <- resR$NPLSDAloadingsperMode3[["item:1"]][[1]]
  .expect_slice_block_ncomp(sw1[["Block.X"]], n, p, q, ncomp)
  .expect_slice_block_ncomp(sw1[["Block.Y"]], n, p, q, ncomp)
  
  # exact unfolded X labels
  exp_unfold_X <- as.vector(outer(dimnames(X)[[2]], dimnames(X)[[3]],
                                  function(v,t) paste0(v, "__t", t)))
  # classical global:
  glX <- resC$NPLSDAloadings$X
  expect_identical(rownames(glX), exp_unfold_X)
  
  # class-2 global:
  exp_unfold_Y <- as.vector(outer(dimnames(Y3)[[2]], dimnames(Y3)[[3]],
                                  function(v,t) paste0(v, "__t", t)))
  glBX <- resR$NPLSDAloadings[["Block.X"]]
  glBY <- resR$NPLSDAloadings[["Block.Y"]]
  expect_identical(rownames(glBX), exp_unfold_X)
  expect_identical(rownames(glBY), exp_unfold_Y)
  
  # slice-wise (first slice)
  expect_identical(rownames(sw1[["Block.X"]]), dimnames(X)[[2]])
  expect_identical(rownames(sw1[["Block.Y"]]), dimnames(Y3)[[2]])
})


