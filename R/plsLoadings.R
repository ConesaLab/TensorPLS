#' Compute NPLS Loadings 
#'
#' Calculates variable loadings using a generalized algorithm for 2–10 components.  
#' Identical math to the original, but without hard-coded steps.
#'
#' @param X 3D array (n × p × k), where n = samples, p = variables, k = time
#' @param Y Matrix (n × q) for classic NPLS-DA, or 3D array (n × q × k) for regression class 2
#' @param ncomp Number of components to extract (2–10, unified algorithm)
#' @param outcome.Y Outcome vector for regression class 2 (required if Y is 3D)
#'
#' @return List of loadings for all requested components
#'
#' 
#'
#' @export
compute_npls_loadings <- function(X, Y, ncomp = 3, outcome.Y = NULL) {
  
  # Validazione input unificata
  .validate_loadings_input_unified(X, Y, ncomp, outcome.Y)
  
  cat("Computing NPLS loadings for ", ncomp, "components...\n")
  
  # checking dim(Y)[2] 
  Ydim <- dim(Y)
  Y_is_3way <- (length(Ydim) == 3L)
  is_regression_class2 <- (length(Ydim) >= 2 && Ydim[2] > 1)
  
  if (is_regression_class2) {
    cat("Performing NPLS-DA regression class 2\n")
  } else {
    cat("Performing Classic NPLS-DA\n")
  }
  
 
  

  ## PARTE 1: ANALISI mixOmics GLOBALE

  
  global_results <- .compute_global_loadings(X, Y, outcome.Y, ncomp, is_regression_class2)
  NPLSDAloadings <- global_results$loadings
  transformedNPLSDAloadings <- global_results$transformed
  

  ## PARTE 2: ANALISI SLICE-WISE per Mode3 

  
  NPLSDAloadingsperMode3 <- .compute_slice_loadings(X, Y, outcome.Y, ncomp, is_regression_class2)
  

  ## PARTE 3: CONSENSUS UNIFICATO 

  
  NPLSDAConsensusloadings <- .compute_consensus_loadings_unified(
    NPLSDAloadingsperMode3, ncomp, is_regression_class2, Y_is_3way, X
  )
  

  ## RETURN RISULTATI

  
  result <- list(
    NPLSDAloadings = NPLSDAloadings,
    transformedNPLSDAloadings = transformedNPLSDAloadings,
    NPLSDAloadingsperMode3 = NPLSDAloadingsperMode3,
    NPLSDAConsensusloadings = NPLSDAConsensusloadings,
    ncomp_used = ncomp
  )
  
  cat("NPLS loadings computed with success \n")
  
  class(result) <- c("npls_loadings_unified", "list")
  return(result)
}



#' Analisi mixOmics Globale per Loadings 
#' @noRd
.compute_global_loadings <- function(X, Y, outcome.Y, ncomp, is_regression_class2) {
  
  # Unfold per mixOmics 
  xdimA3D <- dim(X)
  A <- matrix(X, xdimA3D[1], xdimA3D[2] * xdimA3D[3])
  
  if (is_regression_class2) {
    # NPLS-DA regression class 2 
    ydimB3D <- dim(Y)
    B <- matrix(Y, ydimB3D[1], ydimB3D[2] * ydimB3D[3])
    
    plsdaforplotting <- mixOmics::block.plsda(
      X = list(Block.X = A, Block.Y = B),
      Y = outcome.Y,
      ncomp = ncomp
    )
  } else {
    # Classic NPLS-DA 
    if (length(dim(Y)) == 3) {
      if (!requireNamespace("abind", quietly = TRUE)) {
        stop("Package abind richiesto per questa funzione")
      }
      B3D <- Y
      B3D2 <- abind::abind(B3D, B3D, along = 2)  # Duplica lungo Mode2
      ydims <- dim(B3D)
      Y2 <- matrix(B3D2, ydims[1], ydims[2] * ydims[3] * 2)  
      y_response <- Y[, 1, 1]  
    } else {
      Y2 <- Y  # Se Y è 2D, usalo direttamente
      y_response <- Y[, 1]
    }
    
    plsdaforplotting <- mixOmics::plsda(X = A, Y = y_response, ncomp = ncomp)
  }
  
  # Estrai loadings 
  loadings <- plsdaforplotting$loadings
  
  # Trasformazione per plotting 
  transformed <- NULL
  
  tryCatch({
    # SEMPRE comp=c(1,2) come in codice TEDDY, indipendentemente da ncomp
    tabla <- mixOmics::plotVar(plsdaforplotting, plot = FALSE, comp = c(1, 2))
    
    #  extracts  columns 4,1,2 (names, x, y) like TEDDY 
    # tabla4 <- tabla[,c(4,1,2)] 
    transformed <- tabla[, c(4, 1, 2)]
    
    if (ncol(transformed) == 3) {
      colnames(transformed) <- c("names", "x", "y")
    }
    
  }, error = function(e) {
    cat("Warning: plotVar transformation failed! \n")
    transformed <- NULL
  })
  
  return(list(
    loadings = loadings,
    transformed = transformed
  ))
}

#' Analisi Slice-wise per Loadings 
#' @noRd
.compute_slice_loadings <- function(X, Y, outcome.Y, ncomp, is_regression_class2) {
  
  k <- dim(X)[3]
  NPloadings <- list()
  
  for (i in 1:k) {
    if (is_regression_class2) {
      plsdaforplotting <- mixOmics::block.plsda(
        X = list(Block.X = X[, , i], Block.Y = Y[, , i]),
        Y = outcome.Y,
        ncomp = ncomp
      )
    } else {
      # Classic NPLS-DA: replica la logica NPLSDAmod per slice
      if (length(dim(Y)) == 3) {
        if (!requireNamespace("abind", quietly = TRUE)) {
          stop("Package abind required")
        }
        B3D_slice <- Y[, , i, drop = FALSE]  
        B3D2_slice <- abind::abind(B3D_slice, B3D_slice, along = 2)
        y_response <- Y[, 1, 1]  
      } else {
        y_response <- Y[, 1]  
      }
      
      plsdaforplotting <- mixOmics::plsda(X = X[, , i], Y = y_response, ncomp = ncomp)
    }
    
    NPloadingstita <- plsdaforplotting$loadings
    nameloadingdtita <- paste("item:", i, sep = "")
    loadingstitatemp <- list(NPloadingstita)
    NPloadings[[nameloadingdtita]] <- loadingstitatemp
  }
  
  return(NPloadings)
}

#' Validazione Input 
#' @noRd
.validate_loadings_input_unified <- function(X, Y, ncomp, outcome.Y) {
  
  if (length(dim(X)) != 3) {
    stop("X must be a 3D array (n × p × k)")
  }
  
  if (dim(X)[1] != dim(Y)[1]) {
    stop("X and Y must have the same number of observations")
  }
  
  # Supports 2–10 components with the unified algorithm
  if (ncomp < 2 || ncomp > 10) {
    stop("ncomp must be between 2 and 10 (unified algorithm)")
  }
  
  # Check statistical limits
  n <- dim(X)[1]
  p <- dim(X)[2]
  max_statistical <- min(n - 1, p)
  
  if (ncomp > max_statistical) {
    stop("ncomp (", ncomp, ") exceeds statistical limit (", max_statistical, 
         ") based on n=", n, " subjects and p=", p, " variables")
  }
  
  if (anyNA(X) || anyNA(Y)) {
    stop("X and Y cannot contain NA values")
  }
  
  Y_is_3way <- (length(dim(Y)) == 3L)
  
  if (!is.null(outcome.Y) && !Y_is_3way) {
    stop("outcome.Y was provided but Y is not 3D – inconsistency")
  }
  
  if (Y_is_3way && !is.null(outcome.Y)) {
    cat("Mode: NPLS-DA regression class 2\n")
  } else {
    cat("Mode: Classic NPLS-DA\n")
  }
  
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("The mixOmics package is required for this function")
  }
  
  cat("Validation OK: ncomp=", ncomp, " (unified algorithm)\n")
}



#'
#' 
#' 
#' 
#' @noRd
.compute_consensus_loadings_unified <- function(NPloadings, ncomp, is_regression_class2, Y_is_3way, X) {
  
  if (length(NPloadings) == 0) {
    stop("Nessun loading slice-wise disponibile per il consensus")
  }
  
  n_slices <- length(NPloadings)
  ncomp <- as.integer(ncomp)
  local_env <- new.env()
  for (comp in 1:ncomp) {
    assign(paste0("Xcomp", comp, "loadingsmean"), NULL, envir = local_env)
    assign(paste0("Ycomp", comp, "loadingsmean"), NULL, envir = local_env)
  }

  # RACCOLTA LOADINGS 
  for (i in 1:n_slices) {
    loadings_i <- NPloadings[[i]][[1]]
    
    for (comp in 1:ncomp) {
      
      # Nomi delle variabili
      x_var_name <- paste0("Xcomp", comp, "loadingsmean")
      y_var_name <- paste0("Ycomp", comp, "loadingsmean")
      
      if (is_regression_class2) {
        # *** BLOCK.PLSDA STRUCTURE ***
        
        # X loadings (Block 1)
        if (!is.null(loadings_i[[1]]) && ncol(loadings_i[[1]]) >= comp) {
          x_loading <- as.matrix(loadings_i[[1]][, comp])
          
          current_x <- get(x_var_name, envir = local_env)
          new_x <- cbind(current_x, x_loading)
          assign(x_var_name, new_x, envir = local_env)
        }
        
        # Y loadings (Block 2)  
        if (!is.null(loadings_i[[2]]) && ncol(loadings_i[[2]]) >= comp) {
          y_loading <- as.matrix(loadings_i[[2]][, comp])
          
          current_y <- get(y_var_name, envir = local_env)
          new_y <- cbind(current_y, y_loading)
          assign(y_var_name, new_y, envir = local_env)
        }
        
      } else {
        # *** STANDARD PLSDA STRUCTURE ***
        
        # X loadings  
        if (!is.null(loadings_i$X) && ncol(loadings_i$X) >= comp) {
          x_loading <- as.matrix(loadings_i$X[, comp])
          current_x <- get(x_var_name, envir = local_env)
          new_x <- cbind(current_x, x_loading)
          assign(x_var_name, new_x, envir = local_env)
        }
        
        # Y loadings
        if (!is.null(loadings_i$Y) && ncol(loadings_i$Y) >= comp) {
          y_loading <- as.matrix(loadings_i$Y[, comp])
          
          current_y <- get(y_var_name, envir = local_env)
          new_y <- cbind(current_y, y_loading)
          assign(y_var_name, new_y, envir = local_env)
        }
      }
    }
  }
  # CALCOLO CONSENSUS
  Block.XConsensus <- NULL
  Block.YConsensus <- NULL
  
  for (comp in 1:ncomp) {
    
    x_var_name <- paste0("Xcomp", comp, "loadingsmean")
    y_var_name <- paste0("Ycomp", comp, "loadingsmean")
    
    # *** X CONSENSUS ***
    current_x_matrix <- get(x_var_name, envir = local_env)
    
    if (!is.null(current_x_matrix) && ncol(current_x_matrix) > 0) {
      
      #comp1meanloadings3D <- as.matrix(rowMeans(Xcomp1loadingsmean))
      comp_consensus <- as.matrix(rowMeans(current_x_matrix, na.rm = TRUE))
      colnames(comp_consensus) <- paste0("Comp", comp, "loadingsmean3D")
      
      if (is.null(Block.XConsensus)) {
        Block.XConsensus <- comp_consensus
      } else {
        Block.XConsensus <- cbind(Block.XConsensus, comp_consensus)
      }
      
      cat("X Component", comp, ": consensus from", ncol(current_x_matrix), "slice\n")
    } else {
      cat("Warning: X Component", comp, "not  available \n")
    }
    
    # *** Y CONSENSUS ***
    current_y_matrix <- get(y_var_name, envir = local_env)
    
    if (!is.null(current_y_matrix) && ncol(current_y_matrix) > 0) {
      
      comp_consensus <- as.matrix(rowMeans(current_y_matrix, na.rm = TRUE))
      colnames(comp_consensus) <- paste0("YComp", comp, "loadingsmean3D")
      
      if (is.null(Block.YConsensus)) {
        Block.YConsensus <- comp_consensus
      } else {
        Block.YConsensus <- cbind(Block.YConsensus, comp_consensus)
      }
      
      cat("Y Component", comp, ": consensus from", ncol(current_y_matrix), "slice\n")
    } else {
      cat("Warning: Y Component", comp, "not available\n")
    }
  }
  
  if (!is.null(Block.XConsensus) && !is.null(dimnames(X)[[2]])) {
    rownames(Block.XConsensus) <- dimnames(X)[[2]]
  }
  
  result <- list(
    Block.XConsensus = Block.XConsensus,
    Block.YConsensus = Block.YConsensus
  )
  
  cat("Consensus completed with success! \n")
  
  return(result)
}

#nplsLoadingsOriginal = compute_npls_loadings(X = fullarrayGCTOFX, Y = outcomedummyarray136,ncomp = 3)
#nplsLoadingsOriginal$NPLSDAloadings
#npls
#save(nplsLoadingsOriginal,file = "/Users/alessandrogiordano/Desktop/TEDDY/TestingAndStudying/TestingNPLSDA/NPLSDAPackage/Loadings/Objects/ObjectsUnified/nplsLoadingsOriginalGCTOFX.RData")
