#' N-way Partial Least Squares Discriminant Analysis - 
#'
#' Compute performance metrics for NPLS-DA. Supporta 2-10 components, multivariate Y,
#' and robust error handling
#'
#' @param X Array 3D (n×p×k): eg. subjects, p features, k times
#' @param Y Matrice (n×q) or n x p x q where dim(q) == 1 classic NPLS-DA or  Array 3D (n × p x q) for class2 regression with q > 1
#' @param ncomp Number of components (2-10)
#' @param outcome.Y Outcome vector for  regression class 2 (mandatory if Y is from regression Class2)
#' @param slice_vip Boolean,  about computing nplsda_vips separately for each slice
#'
#'
#'
#' @return List with nplsda_vips, Q², explained variance and predictions for each component
#' 
#' @export
nplsda_vips <- function(X, Y, 
                           ncomp = 3, 
                           outcome.Y = NULL,
                           slice_vip = FALSE) {
  
  # Validation
  .validate(X, Y, ncomp, outcome.Y)
  
  message("Compute NPLS-DA for", ncomp, "components...\n")
  
  # Type of NPLSDA
  Ydim <- dim(Y)
  Y_is_3way <- (length(Ydim) == 3L)
  
  is_regression_class2 <- (Y_is_3way && Ydim[2] > 1)   
  if (is_regression_class2) {
    message("Performing NPLS-DA regression class 2\n")
  } else {
    message("Performing Classic NPLS-DA\n")
  }
  
  
  
  # Setup Sizes
  n <- dim(X)[1]; p <- dim(X)[2]; k <- dim(X)[3]
  Xmat <- matrix(aperm(X, c(1, 2, 3)), n)
  
  # SETUP Y
  
  message("Setup Y \n")
  
  Ymat <- if (Y_is_3way) 
    matrix(aperm(Y, c(1, 2, 3)), n)
  else 
    as.matrix(Y)
  
  message("Sizes of Y:", dim(Ymat), "\n")
  
  if (is_regression_class2) {
    if (is.null(outcome.Y)) {
      stop("outcome.Y required for regression class 2")
    }
    y_response_factor <- factor(outcome.Y)
    message("Regression class 2 detected with", length(levels(y_response_factor)), "groups\n")
  } else {
    # Classic: extraction outcome from Y
    if (Y_is_3way) {
      y_response_factor <- factor(Y[, 1, 1])
    } else {
      y_response_factor <- factor(Y[, 1])
    }
    message("Classic NPLS-DA with", length(levels(y_response_factor)), "groups\n")
  }
  
  #Pls Core  
  message("Execution PLS...\n")
  
  pls_global <- pls_reg(Xmat, Ymat, ncomp = ncomp, cv = TRUE)
  
  message("PLS globale completato senza problemi\n")
  
  # Vip calculations aligned with reference function
  
  VIP2D <- pls_global$VIP
  var_names <- dimnames(X)[[2]]  
  time_names <- dimnames(X)[[3]] 
  if (is.null(var_names)) var_names <- paste0("Var", 1:p)
  if (is.null(time_names)) time_names <- paste0("Time", 1:k)
  vip2d_rownames <- as.vector(outer(var_names, time_names, paste, sep = "_"))
  rownames(VIP2D) <- vip2d_rownames
  colnames(VIP2D) <- paste0("Comp", 1:ncomp)
  VIParray <- array(
    VIP2D,
    dim = c(p, k, ncomp),
    dimnames = list(var_names, time_names, paste0("Comp", 1:ncomp))
  )
  #VIP3DModel1
  if (!slice_vip) {
    
    slice_idx <- seq_len(min(ncomp, k))
    VIP3Dmodel1 <- apply(VIParray[, slice_idx, , drop = FALSE], c(1, 3), mean)
    rownames(VIP3Dmodel1) <- dimnames(X)[[2]]
    colnames(VIP3Dmodel1) <- paste0("t", seq_len(ncomp), "mean3D")
  } 
  else {
    
    VIP_slice <- lapply(seq_len(k), function(tt) {
      X_tt <- matrix(X[, , tt], nrow = n)
      Y_tt <- if (Y_is_3way) matrix(Y[, , tt], nrow = n) else Ymat
      pls_reg(X_tt, Y_tt, ncomp, cv = FALSE)$VIP
    })
    slice_idx <- seq_len(min(ncomp, k))
    VIP3Dmodel1 <- Reduce(`+`, VIP_slice[slice_idx]) / length(slice_idx)
    rownames(VIP3Dmodel1) <- dimnames(X)[[2]]
    colnames(VIP3Dmodel1) <- paste0("t", seq_len(ncomp), "mean3D")
  }
  
  # VIP3DModel2
  if (!slice_vip) {
    
    VIP3Dmodel2 <- apply(VIParray, c(1, 2), sum)
    rownames(VIP3Dmodel2) <- dimnames(X)[[2]]
    colnames(VIP3Dmodel2) <- dimnames(X)[[3]]
  } 
  else {
    # *** slice_vip = TRUE: USA TUTTE LE K SLICE COME NELL'ORIGINALE ***
    VIP3D_list_M2 <- lapply(seq_len(k), function(tt) {
      X_tt <- matrix(X[, , tt], nrow = n)
      Y_tt <- if (Y_is_3way) matrix(Y[, , tt], nrow = n) else Ymat
      rowSums(pls_reg(X_tt, Y_tt, ncomp, cv = FALSE)$VIP)
    })
    VIP3Dmodel2 <- do.call(cbind, VIP3D_list_M2)
    rownames(VIP3Dmodel2) <- dimnames(X)[[2]]
    colnames(VIP3Dmodel2) <- dimnames(X)[[3]]
  }
  
  # Q² 3D CALCULATIONS 
  
  Q2_list <- vector("list", k)
  Q2cum_list <- vector("list", k)
  
  for (tt in seq_len(k)) {
    X_tt <- matrix(X[, , tt], nrow = n)
    
    if (Y_is_3way) {
      Y_tt <- matrix(Y[, , tt], nrow = n)
    } else {
      Y_tt <- Ymat
    }
    
    tryCatch({
      pls_slice <- pls_reg(X_tt, Y_tt, ncomp = ncomp, cv = TRUE)
      
      if (ncol(pls_slice$Q2) == 1) {
        Q2_list[[tt]] <- pls_slice$Q2[, 1]
        Q2cum_list[[tt]] <- pls_slice$Q2cum[, 1]
      } else {
        Q2_list[[tt]] <- rowMeans(pls_slice$Q2)
        Q2cum_list[[tt]] <- rowMeans(pls_slice$Q2cum)
      }
    }, error = function(e) {
      message("Warning: Q² slice", tt, "failed\n")
      Q2_list[[tt]] <<- rep(0, ncomp)
      Q2cum_list[[tt]] <<- rep(0, ncomp)
    })
  }
  
  valid_Q2 <- !sapply(Q2_list, function(x) all(x == 0))
  
  if (sum(valid_Q2) > 0) {
    Q2mean3D_vec <- rowMeans(do.call(cbind, Q2_list[valid_Q2]))
    Q2cummean3D_vec <- rowMeans(do.call(cbind, Q2cum_list[valid_Q2]))
  } else {
    Q2mean3D_vec <- rep(0, ncomp)
    Q2cummean3D_vec <- rep(0, ncomp)
  }
  
  Q2mean3D <- matrix(Q2mean3D_vec, ncol = 1,
                     dimnames = list(paste0("t", seq_len(ncomp)), "Q2mean3D"))
  Q2cummean3D <- matrix(Q2cummean3D_vec, ncol = 1,
                        dimnames = list(paste0("t", seq_len(ncomp)), "Q2cummean3D"))
  
  #   # EXPLAINED VARIANCE 
  
  R2_slices <- vector("list", k)
  
  for (tt in seq_len(k)) {
    tryCatch({
      Xc <- scale(X[, , tt], center = TRUE, scale = TRUE)
      
      if (Y_is_3way) {
        Yc <- scale(Y[, , tt], center = TRUE, scale = TRUE)
      } else {
        Yc <- scale(Ymat, center = TRUE, scale = TRUE)
      }
      
      plsl <- pls_reg(Xc, Yc, ncomp = ncomp, cv = FALSE)
      
      SST_X <- sum(Xc^2)
      SST_Y <- sum(Yc^2)
      
      Xm_hat <- matrix(0, n, p)
      Ym_hat <- matrix(0, n, ncol(Yc))
      SSE_X <- SSE_Y <- numeric(ncomp)
      
      for (h in seq_len(ncomp)) {
        Xm_hat <- Xm_hat + plsl$x.scores[, h] %*% t(plsl$x.loads[, h])
        Ym_hat <- Ym_hat + plsl$x.scores[, h] %*% t(plsl$y.loads[, h])
        SSE_X[h] <- sum((Xc - Xm_hat)^2)
        SSE_Y[h] <- sum((Yc - Ym_hat)^2)
      }
      
      R2X_cum <- 1 - SSE_X / SST_X
      R2Y_cum <- 1 - SSE_Y / SST_Y
      R2X_inc <- c(R2X_cum[1], diff(R2X_cum))   
      R2Y_inc <- c(R2Y_cum[1], diff(R2Y_cum))
      
      R2_slices[[tt]] <- list(R2X = R2X_inc, R2Y = R2Y_inc)
      
    }, error = function(e) {
      message("Warning: Explained variance slice", tt, "failed\n")
      R2_slices[[tt]] <<- list(R2X = rep(0, ncomp), R2Y = rep(0, ncomp))
    })
  }
  
  valid_R2 <- !sapply(R2_slices, function(x) all(x$R2X == 0))
  
  if (sum(valid_R2) > 0) {
    R2X_mat <- do.call(cbind, lapply(R2_slices[valid_R2], `[[`, "R2X"))
    R2Y_mat <- do.call(cbind, lapply(R2_slices[valid_R2], `[[`, "R2Y"))
    
    R2Xmean3D <- rowMeans(R2X_mat)
    R2Ymean3D <- rowMeans(R2Y_mat)
  } else {
    R2Xmean3D <- rep(0, ncomp)
    R2Ymean3D <- rep(0, ncomp)
  }
  
  Explvar3D <- cbind(R2Xmean3D, R2Ymean3D)
  rownames(Explvar3D) <- paste0("t", 1:ncomp)
  colnames(Explvar3D) <- c("R2Xmean3D", "R2Ymean3D")
  
  # EXPLAINED VARIANCE 2D GLOBALE 
  
  Xsc <- scale(Xmat)
  Ysc <- scale(Ymat)
  
  SST_Xg <- sum(Xsc^2)
  SST_Yg <- sum(Ysc^2)
  
  Xm_hat_g <- matrix(0, n, p * k)
  Ym_hat_g <- matrix(0, n, ncol(Ymat))  
  
  R2Xg <- R2Yg <- numeric(ncomp)
  
  for (h in seq_len(ncomp)) {
    Xm_hat_g <- Xm_hat_g + pls_global$x.scores[, h] %*% t(pls_global$x.loads[, h])
    Ym_hat_g <- Ym_hat_g + pls_global$x.scores[, h] %*% t(pls_global$y.loads[, h])
    
    R2Xg[h] <- 1 - sum((Xsc - Xm_hat_g)^2) / SST_Xg
    R2Yg[h] <- 1 - sum((Ysc - Ym_hat_g)^2) / SST_Yg
  }
  
  R2Xg_inc <- c(R2Xg[1], diff(R2Xg))   
  R2Yg_inc <- c(R2Yg[1], diff(R2Yg))
  R2Xg_cum <- R2Xg
  R2Yg_cum <- R2Yg
  
  NPLSDAexplVar <- matrix(
    R2Yg_inc,
    ncol = 1,
    dimnames = list(paste0("t", seq_len(ncomp)), "R2.Y")
  )
  
  ExplVar_global <- cbind(
    R2X = R2Xg_inc,
    R2Xcum = R2Xg_cum,
    R2Y = R2Yg_inc,
    R2Ycum = R2Yg_cum
  )
  rownames(ExplVar_global) <- paste0("t", seq_len(ncomp))
  
  # OUTPUT 
  
  result <- structure(list(
    VIP2D = VIP2D,
    VIP3Dmodel1 = VIP3Dmodel1,
    VIP3Dmodel2 = VIP3Dmodel2,
    Q2 = pls_global$Q2,
    Q2cum = pls_global$Q2cum,
    NPLSDAQ2mean3D = Q2mean3D,
    NPLSDAQ2cummean3D = Q2cummean3D,
    y.pred = pls_global$y.pred,
    Explvar3D = Explvar3D,
    NPLSDAexplVar = NPLSDAexplVar,
    explvar = ExplVar_global,
    residuals = pls_global$resid,
    ncomp_used = ncomp,
    regression_class2 = is_regression_class2
  ), class = "nplsda_vips")
  
  message("NPLS-DA with nplsda_vips completed with Success!\n")
  
  return(result)
}

#' Validate Inputs for Unified NPLS-DA
#'
#' This function validates the input arguments for the unified NPLS-DA function. 
#'
#' @param X A 3D array of size \eqn{n × p × k}, where \eqn{n} is the number of subjects,
#'   \eqn{p} is the number of variables, and \eqn{k} is an additional mode (e.g., time or block).
#' @param Y A matrix (\eqn{n × q}) or a 3D array (\eqn{n × q × k}) containing the response(s).
#' @param ncomp Integer. The number of components to extract. Must be between 2 and 10, 
#'   and also respect the statistical limit \eqn{\min(n-1, p)}.
#' @param outcome.Y Optional. An outcome variable for regression/class 2 mode. 
#'   If provided, Y must be a 3D array.
#'
#'
#' @return This function does not return a value. It stops execution if validation fails, 
#' or prints messages confirming successful validation and the analysis mode.
#'
#' @details
#' The function checks:
#' \itemize{
#'   \item X must be a 3D array.
#'   \item X and Y must have the same number of subjects.
#'   \item ncomp must be between 2 and 10 and not exceed \eqn{\min(n-1, p)}.
#'   \item Neither X nor Y may contain NA values.
#'   \item If \code{outcome.Y} is provided, Y must be a 3D array (regression class 2 mode).
#' }
#'
#' Requires the \pkg{mixOmics} package.
#'
#' @noRd
.validate <- function(X, Y, ncomp, outcome.Y) {
  
  # Reuse validation identical to loadings
  if (length(dim(X)) != 3) {
    stop("X should be a 3D array (n × p × k)")
  }
  
  if (dim(X)[1] != dim(Y)[1]) {
    stop("X and Y must have the same number of observations")
  }
  
  if (ncomp < 2 || ncomp > 10) {
    stop("ncomp must be between 2 and 10 (unified algorithm)")
  }
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  max_statistical <- min(n - 1, p)
  
  if (ncomp > max_statistical) {
    stop("ncomp (", ncomp, 
         ") exceeds the statistical limit (", max_statistical, 
         ") based on n = ", n, " subjects and p = ", p, " variables")
  }
  
  
  if (anyNA(X) || anyNA(Y)) {
    stop("X and Y cannot contain NA values")
  }
  
  # Checking outcome.Y consistency
  Y_is_3way <- (length(dim(Y)) == 3L)
  
  if (!is.null(outcome.Y) && !Y_is_3way) {
    stop("outcome.Y was provided but Y is not 3D - inconsistency detected")
  }
  
  if (Y_is_3way && !is.null(outcome.Y)) {
    message("Mode: NPLS-DA regression class 2 - performance metrics\n")
  } else {
    message("Mode: Classic NPLS-DA - performance metrics\n")
  }
  
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package 'mixOmics' is required for this function")
  }
  
  message("Validation OK: ncomp = ", ncomp, " (unified performance metrics)\n")
}

