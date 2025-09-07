# tests/testthat/test_compute_npls_variates.R

library(testthat)

skip_on_cran()
skip_if_not_installed("mixOmics")

set.seed(123)

# ---- Small helpers ----------------------------------------------------------

# Generate 3D predictor tensor X with a class effect on a subset of variables/times
gen_X <- function(n = 30, p = 6, k = 3, effect = 1.2, groups = NULL) {
  X <- array(rnorm(n * p * k, 0, 1), dim = c(n, p, k))
  if (!is.null(groups)) {
    idx <- groups == levels(groups)[1]  # push class 1 positively on first half of features
    p_sig <- max(2, floor(p / 2))
    k_sig <- max(1, floor(k / 2))
    for (tt in seq_len(k_sig)) {
      X[idx, 1:p_sig, tt]  <- X[idx, 1:p_sig, tt]  + effect
      X[!idx, 1:p_sig, tt] <- X[!idx, 1:p_sig, tt] - effect
    }
  }
  X
}

# Classic-mode Y: matrix n x q (q>=1). First column encodes outcome (numeric 1/2 is fine;
# your function should convert to factor internally).
gen_Y_classic <- function(groups, q = 1) {
  n <- length(groups)
  Y <- matrix(0, nrow = n, ncol = q)
  Y[, 1] <- as.integer(groups) # 1/2
  if (q > 1) Y[, 2:q] <- matrix(rnorm(n * (q - 1), 0, 1), n, q - 1)
  Y
}

# Class-2 Y: 3D array n x q x k (q>1) with mild correlation to groups
gen_Y_class2 <- function(n = 30, q = 2, k = 3, groups = NULL, effect = 0.8) {
  Y <- array(rnorm(n * q * k, 0, 1), dim = c(n, q, k))
  if (!is.null(groups)) {
    idx <- groups == levels(groups)[1]
    for (tt in seq_len(k)) {
      Y[idx,  , tt] <- Y[idx,  , tt] + effect
      Y[!idx, , tt] <- Y[!idx, , tt] - effect
    }
  }
  Y
}
.only_block_variates <- function(vars) {
  # Se è una lista di variates (block.plsda), rimuovi l'eventuale componente "Y".
  if (is.list(vars)) {
    nm <- names(vars)
    if (!is.null(nm)) {
      vars <- vars[nm != "Y"]
    }
  }
  vars
}


# ---- Input validation -------------------------------------------------------

test_that("Input validation: shape, ncomp, NA, and outcome consistency", {
  n <- 20; p <- 4; k <- 2
  groups <- factor(sample(c("A", "B"), n, TRUE))
  X <- gen_X(n, p, k, groups = groups)
  Y_classic <- gen_Y_classic(groups, q = 1)
  
  # X must be 3D
  expect_error(
    compute_npls_variates(X = matrix(0, n, p), Y = Y_classic, ncomp = 2),
    "X should be an array 3D (n × p × k)", fixed = TRUE
  )
  
  # X and Y same n
  expect_error(
    compute_npls_variates(X = X, Y = Y_classic[-1, , drop = FALSE], ncomp = 2),
    "same numner of observations", fixed = TRUE
  )
  
  # ncomp >= 2
  expect_error(
    compute_npls_variates(X = X, Y = Y_classic, ncomp = 1),
    "ncomp should be >=2", fixed = TRUE
  )
  
  # NA not allowed
  Xna <- X; Xna[1,1,1] <- NA
  expect_error(
    compute_npls_variates(X = Xna, Y = Y_classic, ncomp = 2),
    "Needs to do Imputation", fixed = TRUE
  )
  
  # outcome.Y provided but Y is NOT 3D -> error
  expect_error(
    compute_npls_variates(X = X, Y = Y_classic, ncomp = 2, outcome.Y = groups),
    "outcome.Y given but Y isn't 3D", fixed = TRUE
  )
})

# ---- Classic mode (Y as matrix) --------------------------------------------

test_that("Classic mode: returns coherent structures, shapes, and consensus", {
  n <- 28; p <- 6; k <- 3; ncomp <- 3
  groups <- factor(sample(c("A", "B"), n, TRUE))
  X <- gen_X(n, p, k, effect = 1.3, groups = groups)
  Y <- gen_Y_classic(groups, q = 2)   # q>=1 ok; first column is used as outcome
  
  res <- compute_npls_variates(X, Y, ncomp = ncomp)
  
  # Result container
  expect_s3_class(res, "npls_variates")
  expect_true(all(c("NPLSDAvariates", "NPLSDAvariatesperMode3", "NPLSDAConsensusvariates") %in% names(res)))
  
  # Global variates (mixOmics::plsda): list with $X and $Y matrices
  expect_true(is.list(res$NPLSDAvariates))
  expect_true(all(c("X", "Y") %in% names(res$NPLSDAvariates)))
  expect_equal(nrow(res$NPLSDAvariates$X), n)
  expect_gte(ncol(res$NPLSDAvariates$X), ncomp)
  expect_equal(nrow(res$NPLSDAvariates$Y), n)
  
  # Slice-wise list length == k, each entry is list(variates)
  expect_length(res$NPLSDAvariatesperMode3, k)
  expect_true(all(grepl("^item:", names(res$NPLSDAvariatesperMode3))))
  v1 <- res$NPLSDAvariatesperMode3[[1]][[1]]
  expect_true(all(c("X", "Y") %in% names(v1)))
  expect_equal(nrow(v1$X), n)
  expect_gte(ncol(v1$X), ncomp)
  
  # Consensus has up to ncomp columns, n rows
  cons <- res$NPLSDAConsensusvariates
  expect_true(all(c("Block.XConsensus", "Block.YConsensus") %in% names(cons)))
  expect_equal(nrow(cons$Block.XConsensus), n)
  expect_lte(ncol(cons$Block.XConsensus), ncomp)
  expect_equal(nrow(cons$Block.YConsensus), n)
  expect_lte(ncol(cons$Block.YConsensus), ncomp)
})

# ---- Regression class-2 mode (Y as 3D with q>1, outcome.Y provided) --------

test_that("Regression class-2 mode: runs and returns block variates + consensus", {
  n <- 30; p <- 5; k <- 3; q <- 2; ncomp <- 3
  groups <- factor(sample(c("A", "B"), n, TRUE))
  X <- gen_X(n, p, k, effect = 1.1, groups = groups)
  Y3 <- gen_Y_class2(n = n, q = q, k = k, groups = groups, effect = 0.7)
  
  res <- compute_npls_variates(X, Y3, ncomp = ncomp, outcome.Y = groups)
  
  # ---- Global block variates (filtrati)
  vv <- .only_block_variates(res$NPLSDAvariates)
  expect_true(is.list(vv))
  expect_equal(length(vv), 2L)                      # ora passa
  expect_true(all(vapply(vv, is.matrix, logical(1))))
  expect_true(all(vapply(vv, nrow, integer(1)) == n))
  expect_true(all(vapply(vv, ncol, integer(1)) >= ncomp))
  
  # ---- Slice-wise (filtra anche qui)
  v1 <- .only_block_variates(res$NPLSDAvariatesperMode3[[1]][[1]])
  expect_true(is.list(v1) && length(v1) == 2L)
  expect_true(all(vapply(v1, is.matrix, logical(1))))
  expect_true(all(vapply(v1, nrow, integer(1)) == n))
  
  # ---- Consensus: invariato
  cons <- res$NPLSDAConsensusvariates
  expect_true(all(c("Block.XConsensus", "Block.YConsensus") %in% names(cons)))
  expect_equal(nrow(cons$Block.XConsensus), n)
  expect_lte(ncol(cons$Block.XConsensus), ncomp)
  expect_equal(nrow(cons$Block.YConsensus), n)
  expect_lte(ncol(cons$Block.YConsensus), ncomp)
})


# ---- Determinism 

test_that("Deterministic results given same seed and data", {
  n <- 24; p <- 5; k <- 3; ncomp <- 2
  groups <- factor(sample(c("A", "B"), n, TRUE))
  X <- gen_X(n, p, k, effect = 1.0, groups = groups)
  Y <- gen_Y_classic(groups, q = 1)
  
  set.seed(777)
  r1 <- compute_npls_variates(X, Y, ncomp = ncomp)
  set.seed(777)
  r2 <- compute_npls_variates(X, Y, ncomp = ncomp)
  
  # Compare consensus X scores; allow tiny numerical noise
  expect_equal(r1$NPLSDAConsensusvariates$Block.XConsensus,
               r2$NPLSDAConsensusvariates$Block.XConsensus,
               tolerance = 1e-8)
})

