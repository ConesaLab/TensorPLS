#' Compute NPLS Variates
#' 
#'
#' @param X Array 3D (n × p × k) dove n = subject, p = variables, k = time
#' @param Y Matrix (n × q) for classic NPLS-DA or 3D Array (n × q × k) for class 2 regression
#' @param ncomp Number of components to extract (2, 3, o 4)
#' @param outcome.Y Outcome vector for regression class 2 
#' 
#' @return List with:
#'   \item{NPLSDAvariates}{Variates for individual global plot}
#'   \item{NPLSDAvariatesperMode3}{Variates for each temporal slice}
#'   \item{NPLSDAConsensusvariates}{Consensus variates}
#' 
#'
#' 
#' @note
#' For regression class 2 (Y 3D), the parameter outcome.Y è obbligatorio.
#' 
#' @examples
#' \dontrun{
#' # Classic NPLS-DA
#' variates <- compute_npl_variates_gen(
#'   X = my_3d_array,
#'   Y = my_outcome_matrix,
#'   ncomp = 3
#' )
#' 
#' # Plot individuali
#' mixOmics::plotIndiv(variates$NPLSDAvariates)
#' 
#' # Regression class 2
#' variates <- compute_npl_variates_gen(
#'   X = my_3d_array,
#'   Y = my_3d_response,
#'   outcome.Y = my_groups,
#'   ncomp = 3
#' )
#' }
#' 
#' @seealso \code{\link{compute_npls_factors}}, \code{\link{mixOmics::plotIndiv}}
#' 
#' @export
compute_npls_variates <- function(X, Y, ncomp = 3, outcome.Y = NULL) {
  
  # Validazione input
  .validate_variates_input(X, Y, ncomp, outcome.Y)
  
  message("Computing NPLS variates...\n")
  
  Ydim <- dim(Y)
  Y_is_3way <- (length(Ydim) == 3L)
  is_regression_class2 <- (Y_is_3way && Ydim[2] > 1)   
  
  if (is_regression_class2) {
    message("Performing NPLS-DA regression class 2\n")
  } else {
    message("Performing Classic NPLS-DA\n")
  }
  
  NPLSDAvariates <- .compute_global_variates(X, Y, outcome.Y, ncomp, is_regression_class2)

  #  SLICE-WISE for Mode3

  
  NPLSDAvariatesperMode3 <- .compute_slice_variates(X, Y, outcome.Y, ncomp, is_regression_class2)
  
  ## PART 3: CONSENSUS CONFIGURATION 
  
  NPLSDAConsensusvariates <- .compute_consensus_variates_exact(
    NPLSDAvariatesperMode3, ncomp, is_regression_class2, Y_is_3way, X
  )
  

  ## RETURN result 

  
  result <- list(
    NPLSDAvariates = NPLSDAvariates,
    NPLSDAvariatesperMode3 = NPLSDAvariatesperMode3,
    NPLSDAConsensusvariates = NPLSDAConsensusvariates
  )
  
  message("NPLS variates computed with success!\n")
  
  class(result) <- c("npls_variates", "list")
  return(result)
}

#' Validation Input for Variates 
#' @noRd
.validate_variates_input <- function(X, Y, ncomp, outcome.Y) {
  
  if (length(dim(X)) != 3) {
    stop("X should be an array 3D (n × p × k)")
  }
  
  if (dim(X)[1] != dim(Y)[1]) {
    stop("X and Y needs to have the same numner of observations")
  }
  
  if (ncomp < 2) {
    stop("ncomp should be >=2")
  }
  
  
  if (anyNA(X) || anyNA(Y)) {
    stop("Needs to do Imputation")
  }
  
  Y_is_3way <- (length(dim(Y)) == 3L)
  
  # Solo se outcome.Y è fornito esplicitamente, allora è regression class 2
  if (!is.null(outcome.Y) && !Y_is_3way) {
    stop("outcome.Y given but Y isn't 3D")
  }
  
  #classic NPLS-DA 
  if (Y_is_3way && !is.null(outcome.Y)) {
    message("Modality: NPLS-DA regression class 2\n")
  } else {
    message("Modality: Classic NPLS-DA\n")
  }
  
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package mixOmics required")
  }
}

#' Mixomics
#' @noRd
.compute_global_variates <- function(X, Y, outcome.Y, ncomp, is_regression_class2) {
  
  
  xdimA3D <- dim(X)
  A <- matrix(X, xdimA3D[1], xdimA3D[2] * xdimA3D[3])
  n <- nrow(A)
  if (is_regression_class2) {
    # NPLS-DA regression class 2
    ydimB3D <- dim(Y)
    B <- matrix(Y, ydimB3D[1], ydimB3D[2] * ydimB3D[3])
    ncomp_eff <- min(ncomp, n - 1, ncol(A), ncol(B))
    plsdaforplotting <- mixOmics::block.plsda(
      X = list(Block.X = A, Block.Y = B),
      Y = outcome.Y,
      ncomp = ncomp_eff
    )
  } else {
    ncomp_eff <- min(ncomp, n - 1, ncol(A))
    # Classic NPLS-DA 
    if (length(dim(Y)) == 3) {
      if (!requireNamespace("abind", quietly = TRUE)) {
        stop("Package abind required for this function")
      }
      B3D <- Y
      B3D2 <- abind::abind(B3D, B3D, along = 2)  
      ydims <- dim(B3D)
      Y2 <- matrix(B3D2, ydims[1], ydims[2] * ydims[3] * 2)  
      y_response <- Y[, 1, 1]  
    } else {
      Y2 <- Y  
      y_response <- Y[, 1]
    }
    
    plsdaforplotting <- mixOmics::plsda(X = A, Y = y_response, ncomp = ncomp_eff)
  }
  
  
  return(plsdaforplotting$variates)
}

#' Analisi Slice-wise per Variates
#' @noRd
.compute_slice_variates <- function(X, Y, outcome.Y, ncomp, is_regression_class2) {
  
  k <- dim(X)[3]
  NPvariates <- list()
 
  
  for (i in 1:k) {
    Xi <- X[, , i, drop = TRUE]
    n  <- nrow(Xi)
    if (is_regression_class2) {
      n  <- nrow(Xi)
      Yi <- Y[, , i, drop = FALSE][,,1]  # n x q
      ncomp_eff <- min(ncomp, n - 1, ncol(Xi), ncol(Yi))
      plsdaforplotting <- mixOmics::block.plsda(
        X = list(Block.X = Xi, Block.Y = Yi),   
        Y = outcome.Y,
        ncomp = ncomp_eff
      )
    } else {
      Xi <- X[, , i, drop = FALSE][,,1]  # n x p
      ncomp_eff <- min(ncomp, n - 1, ncol(Xi))
      
      if (length(dim(Y)) == 3) {
        if (!requireNamespace("abind", quietly = TRUE)) {
          stop("Package abind required for this function")
        }
        B3D_slice <- Y[, , i, drop = FALSE]  # Mantieni 3D
        B3D2_slice <- abind::abind(B3D_slice, B3D_slice, along = 2)
        Y2_slice <- matrix(B3D2_slice, dim(B3D_slice)[1], dim(B3D_slice)[2] * dim(B3D_slice)[3] * 2)
        y_response <- Y[, 1, 1]  
      } else {
        Y2_slice <- Y  
        y_response <- Y[, 1]
      }
      
      plsdaforplotting <- mixOmics::plsda(X = Xi, Y = y_response, ncomp = ncomp_eff)
    }
    
    NPvariatestita <- plsdaforplotting$variates
    namevariatestita <- paste("item:", i, sep = "")
    variatestitatemp <- list(NPvariatestita)
    NPvariates[[namevariatestita]] <- variatestitatemp
  }
  
  return(NPvariates)
}

.compute_consensus_variates_exact <- function(NPvariates, ncomp, is_regression_class2, Y_is_3way, X) {
  k <- dim(X)[3]  
  
  X_means_cols <- list()
  Y_means_cols <- list()
  
  for (h in seq_len(ncomp)) {
    Xcols <- list()
    Ycols <- list()
    
    for (i in seq_len(k)) {
      variates_i <- NPvariates[[i]][[1]]
      
      if (is_regression_class2) {
        # block.plsda
        if (!is.null(variates_i[[1]]) && ncol(variates_i[[1]]) >= h) {
          Xcols[[length(Xcols) + 1]] <- variates_i[[1]][, h, drop = FALSE]
        }
        if (!is.null(variates_i[[2]]) && ncol(variates_i[[2]]) >= h) {
          Ycols[[length(Ycols) + 1]] <- variates_i[[2]][, h, drop = FALSE]
        }
      } else {
        if (!is.null(variates_i$X) && ncol(variates_i$X) >= h) {
          Xcols[[length(Xcols) + 1]] <- variates_i$X[, h, drop = FALSE]
        }
        if (!is.null(variates_i$Y) && ncol(variates_i$Y) >= h) {
          Ycols[[length(Ycols) + 1]] <- variates_i$Y[, h, drop = FALSE]
        }
      }
    }
    
    if (length(Xcols) > 0) {
      Xmean_h <- matrix(rowMeans(do.call(cbind, Xcols)), ncol = 1)
      colnames(Xmean_h) <- paste0("Comp", h, "mean3D")
      X_means_cols[[length(X_means_cols) + 1]] <- Xmean_h
    }
    if (length(Ycols) > 0) {
      Ymean_h <- matrix(rowMeans(do.call(cbind, Ycols)), ncol = 1)
      colnames(Ymean_h) <- paste0("Comp", h, "mean3D")
      Y_means_cols[[length(Y_means_cols) + 1]] <- Ymean_h
    }
  }
  
  Block.XConsensus <- if (length(X_means_cols) > 0) do.call(cbind, X_means_cols) else NULL
  Block.YConsensus <- if (length(Y_means_cols) > 0) do.call(cbind, Y_means_cols) else NULL
  
  list(
    Block.XConsensus = Block.XConsensus,
    Block.YConsensus = Block.YConsensus
  )
}


