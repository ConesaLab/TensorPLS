#' Compute NPLS Factors with Core Arrays 
#'
#' Computes NPLS factors and core arrays \code{G} and \code{Gu} following
#' the original NPLSDAmod algorithm, with handling of the classical
#' and regression class-2 modes and consistent labeling.
#'
#' @param X A 3D array (n × p × k) where \code{n} = samples, \code{p} = variables,
#'   \code{k} = time/slices.
#' @param Y A response matrix (n × q) for classical NPLS-DA, or a 3D array
#'   (n × q × k) for regression class-2.
#' @param ncomp Integer. Number of components to extract (2–10).
#' @param tol Convergence tolerance for the inner iterative loop. Default: 1e-16.
#' @param max_iter Maximum number of iterations for convergence. Default: 10000.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{FactorsX}: list with Mode1 (scores \code{Tt}), Mode2 (variable loadings \code{WsupraJ}),
#'         Mode3 (time loadings \code{WsupraK})
#'   \item \code{FactorsY}: list with Y factors if \code{Y} is 3D (Mode1 \code{U}, Mode2 \code{QsupraJ},
#'         Mode3 \code{QsupraK}); \code{NULL} otherwise
#'   \item \code{B}: regression coefficient matrix
#'   \item \code{G}: core array (3D, \code{ncomp × ncomp × ncomp})
#'   \item \code{Gu}: list of per-component core matrices (as in the original)
#' }
#'
#' @export
compute_npls_factors <- function(X, Y, ncomp,
                                 tol = 1e-16,
                                 max_iter = 10000) {
  
  .validate_factors_input_unified(X, Y, ncomp)
  
  cat("Computing NPLS factors with core arrays ...\n")
  # Core computation ( identical to the original implementation)
  npls_result <- .npls_factors_with_core_arrays(X, Y, ncomp, tol, max_iter)
  
  Y_is_3way <- (length(dim(Y)) == 3L)
  
  FactorsX <- list(
    Mode1 = npls_result$Tt,
    Mode2 = npls_result$WsupraJ,
    Mode3 = npls_result$WsupraK
  )
  
  FactorsY <- NULL
  if (Y_is_3way) {
    FactorsY <- list(
      Mode1 = npls_result$U,
      Mode2 = npls_result$QsupraJ,
      Mode3 = npls_result$QsupraK
    )
  }
  
  cat("NPLS factors and core arrays computed successfully.\n")
  
  return(list(
    FactorsX = FactorsX,
    FactorsY = FactorsY,
    B = npls_result$B,
    G = npls_result$G,
    Gu = npls_result$Gu
  ))
}

#' Core NPLS algorithm with G and Gu (internal)
#'
#' Implements the original NPLS core loop and computes \code{G} and \code{Gu}
#' during the main iterations.  is identical to the existing code.
#'
#' @noRd
.npls_factors_with_core_arrays <- function(X, Y, ncomp, tol = 1e-16, max_iter = 10000) {
  
  n <- dim(X)[1]; p <- dim(X)[2]; k <- dim(X)[3]
  Y_is_3way <- (length(dim(Y)) == 3L)
  q <- if (Y_is_3way) dim(Y)[2] else ncol(Y)
  
  # Initialization (kept identical)
  Tt <- U <- WsupraJ <- WsupraK <- QsupraJ <- QsupraK <- NULL
  B  <- G <- matrix(0, ncol = ncomp, nrow = ncomp)
  Gu <- vector("list", ncomp)
  
  # Unfolding (kept identical)
  Xmat <- matrix(aperm(X, c(1, 2, 3)), n, p * k)
  if (Y_is_3way) {
    Ymat <- matrix(aperm(Y, c(1, 2, 3)), n, q * k)
  } else {
    # Classical NPLS-DA: replicate Y across slices (as in the original)
    Y_ext <- array(0, dim = c(n, q, k))
    for (i in 1:k) Y_ext[ , , i] <- as.matrix(Y)
    Ymat <- matrix(aperm(Y_ext, c(1, 2, 3)), n, q * k)
  }
  
  # Working copies (X is not deflated; only Y is deflated — identical behavior)
  X_work <- Xmat
  Y_work <- Ymat
  
  # Main component loop (identical )
  for (f in 1:ncomp) {
    
    # Initialize u
    if (ncol(Y_work) > 0) {
      Uf <- svd(Y_work)$u
      u  <- Uf[, 1]
    } else {
      u <- rnorm(n)
    }
    
    # Fixed-point iterations
    for (it in 1:max_iter) {
      
      # 1) Z and SVD
      Zrow   <- t(X_work) %*% u
      Z      <- matrix(Zrow, nrow = p, ncol = k)
      svd_z  <- svd(Z)
      wsupraj <- svd_z$u[, 1]
      wsuprak <- svd_z$v[, 1]
      
      # 2) X scores
      tf <- X_work %*% kronecker(wsuprak, wsupraj)
      
      # 3) V and SVD for Y
      Vrow  <- t(Y_work) %*% tf
      V     <- matrix(Vrow, nrow = q, ncol = k)
      svd_v <- svd(V)
      qsupraj <- svd_v$u[, 1]
      qsuprak <- svd_v$v[, 1]
      
      # 4) Y scores
      uf <- Y_work %*% kronecker(qsuprak, qsupraj)
      
      # Convergence check
      if (sum((uf - u)^2) < tol) {
        cat("Component", f, "converged in", it, "iterations\n")
        break
      }
      u <- uf
    }
    
    # Save factors
    Tt       <- cbind(Tt, tf)
    WsupraJ  <- cbind(WsupraJ, wsupraj)
    WsupraK  <- cbind(WsupraK, wsuprak)
    U        <- cbind(U, uf)
    QsupraJ  <- cbind(QsupraJ, qsupraj)
    QsupraK  <- cbind(QsupraK, qsuprak)
    
    # Regression coefficient (identical)
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("Package MASS is required for ginv().")
    }
    bf <- MASS::ginv(t(Tt) %*% Tt) %*% t(Tt) %*% uf
    B[1:length(bf), f] <- bf
    
    # --- Compute Gu for this component (identical to original) ---
    TM  <- MASS::ginv(t(Tt) %*% Tt) %*% t(Tt)
    WkM <- MASS::ginv(t(WsupraK) %*% WsupraK) %*% t(WsupraK)
    WjM <- MASS::ginv(t(WsupraJ) %*% WsupraJ) %*% t(WsupraJ)
    
    Gu[[f]] <- TM %*% X_work %*% kronecker(t(WkM), t(WjM))
    
    # Deflate Y only
    Y_work <- Y_work - Tt %*% bf %*% t(kronecker(qsuprak, qsupraj))
    
    # Re-init u for next component
    if (f < ncomp) {
      Uf <- svd(Y_work)$u
      u  <- Uf[, 1]
    }
  }
  
  # Build G (3D core) exactly as in the original 
    if (ncomp >= 2) {
    G <- array(as.vector(Gu[[2]]), c(ncomp, ncomp, ncomp))
  } else {
    G <- array(as.vector(Gu[[1]]), c(ncomp, ncomp, ncomp))
  }
  
  # --- Labeling (adds names only; numeric results are unchanged) ---
  comp_names <- paste0("Comp", seq_len(ncomp))
  
  rownames(Tt)      <- dimnames(X)[[1]]
  colnames(Tt)      <- comp_names
  
  rownames(WsupraJ) <- dimnames(X)[[2]]
  colnames(WsupraJ) <- comp_names
  
  rownames(WsupraK) <- dimnames(X)[[3]]
  colnames(WsupraK) <- comp_names
  
  rownames(U)       <- dimnames(X)[[1]]
  colnames(U)       <- comp_names
  
  if (Y_is_3way) {
    rownames(QsupraJ) <- dimnames(Y)[[2]]
    rownames(QsupraK) <- dimnames(Y)[[3]]
  } else {
    rownames(QsupraJ) <- colnames(Y)
    rownames(QsupraK) <- dimnames(X)[[3]]
  }
  colnames(QsupraJ) <- comp_names
  colnames(QsupraK) <- comp_names
  
  rownames(B) <- comp_names
  colnames(B) <- comp_names
  
  dimnames(G) <- list(comp_names, comp_names, comp_names)
  
  return(list(
    WsupraJ = WsupraJ, WsupraK = WsupraK,
    QsupraJ = QsupraJ, QsupraK = QsupraK,
    Tt = Tt, U = U, B = B,
    G = G, Gu = Gu
  ))
}

#'  input validation for factor computation (internal)
#' @noRd
.validate_factors_input_unified <- function(X, Y, ncomp) {
  
  if (length(dim(X)) != 3L) {
    stop("X must be a 3D array (n × p × k).")
  }
  if (dim(X)[1] != dim(Y)[1]) {
    stop("X and Y must have the same number of observations (n).")
  }
  
  # 2–10 components as in the  approach
  if (ncomp < 2L || ncomp > 10L) {
    stop("ncomp must be between 2 and 10")
  }
  
  # Statistical feasibility identical to the loadings validator
  n <- dim(X)[1]; p <- dim(X)[2]
  max_statistical <- min(n - 1L, p)
  if (ncomp > max_statistical) {
    stop("ncomp (", ncomp,
         ") exceeds the statistical limit (", max_statistical,
         ") based on n = ", n, " samples and p = ", p, " variables.")
  }
  
  if (anyNA(X) || anyNA(Y)) {
    stop("X and Y must not contain NA values.")
  }
  
  # Informative message: classical vs regression class-2
  Ydim <- dim(Y)
  Y_is_3way <- (length(Ydim) == 3L)
  if (Y_is_3way && Ydim[2] > 1L) {
    cat("Mode: Regression class-2 (3D Y with q > 1)\n")
  } else {
    cat("Mode: Classical NPLS-DA\n")
  }
  
  cat("Validation OK: ncomp =", ncomp,"\n")
}
#plsFactorsnew = compute_npls_factors(X = fullarrayGeneExpression, Y  = outcomedummyarray136,ncomp = 3)

#cor(plsFactorsnew$FactorsX$Mode1,NPLSDAGeneExpressionFactors$FactorsX$Mode1)
#cor(plsFactorsnew$FactorsX$Mode2,NPLSDAGeneExpressionFactors$FactorsX$Mode2)
#cor(plsFactorsnew$FactorsX$Mode3,NPLSDAGeneExpressionFactors$FactorsX$Mode3)
#cor(plsFactorsnew$FactorsY$Mode1,NPLSDAGeneExpressionFactors$FactorsY$Mode1)
#cor(plsFactorsnew$G,NPLSDAGeneExpressionFactors$G)
#plsFactorsnew$G
#NPLSDAGeneExpressionFactors$G
#plsFactorsnew$B
#NPLSDAGeneExpressionFactors$B
