#' PLS REGRESSION with incremental Q^2 
#' PLS Regression with Incremental cross validation 
#' @param X Matrice predictors (n x p)
#' @param Y Matrice response (n x q)  
#' @param ncomp Numero di componenti PLS (default: 2)
#' @param cv Logico, se TRUE esegue cross-validation (default: TRUE)
#'
#' @details
#' Q² Computing:
#' Q²[h] = 1 - (PRESS[h] / RSS[h])
#'  RSS[h] are residuals after h-1 components 
#'
#' Q²CUM : Q²cum[h] = 1 - ∏(1 - Q²[j]) per j=1:h
#'
#' @return list with components PLS, metrics of validations and predictions
pls_reg <- function(X, Y, ncomp = 2, cv = TRUE) {
  
  # HELPER functions
  .check_matrix <- function(A, name) {
    if (!is.matrix(A)) stop(sprintf("%s should be a numeric matrix", name))
    if (anyNA(A))      stop(sprintf("%s contain NA – needs imputation", name))
    storage.mode(A) <- "double"
  }
  
  .kfold <- function(n, k = 10) {
    id <- rep_len(1:k, n)[sample.int(n)]
    lapply(1:k, function(i) which(id == i))
  }
  
  .vip <- function(W, R2y) {
    p <- nrow(W); h <- ncol(W)
    Rd <- matrix(0, h, h)
    for (j in 1:h) Rd[1:j, j] <- R2y[1:j]
    sqrt((W^2) %*% Rd %*% diag(p/cumsum(R2y), h))
  }
  
  # Input validation
  .check_matrix(X, "X"); .check_matrix(Y, "Y")
  if (nrow(X) != nrow(Y)) stop("X and Y: #differents rows")
  if (ncomp < 2)          stop("ncomp should be ≥2")
  if (!is.logical(cv))    cv <- FALSE
  
  n  <- nrow(X); p <- ncol(X); q <- ncol(Y)
  nc <- min(ncomp, n - 1, p)
  
  # Preprocessing
  Xs <- scale(X); Ys <- scale(Y)
  
  # Component matrix initialization
  W  <- matrix(0,  p, nc)    # Raw weights
  Tt <- matrix(0,  n, nc)    # X-scores  
  U  <- matrix(0,  n, nc)    # Y-scores
  Cc <- matrix(0,  q, nc)    # Y-loadings
  Pp <- matrix(0,  p, nc)    # X-loadings
  
  # Inizialization CROSS-VALIDATION incremental
  if (cv) {
    folds  <- .kfold(n, 10)
    RSS    <- matrix(0, nc + 1, q)    
    RSS[1, ] <- rep(n - 1, q)         
    PRESS  <- matrix(0, nc, q)         
    Q2     <- matrix(0, nc, q)        
  }
  
  # Copies for deflation 
  Xs_deflated <- Xs
  Ys_deflated <- Ys
  
  # Algorithm PLS with incremental CROSS-VALIDATION 
  for (h in 1:nc) {
    
    # Pls algorithm standard (NIPALS)
    u <- Ys_deflated[, 1]; w <- rep(1, p)
    
    for (iter in 1:100) {
      # Step 1: Weight vector from Y-scores
      w_new <- t(Xs_deflated) %*% u / sum(u^2)
      w_new <- w_new / sqrt(sum(w_new^2))
      
      # Step 2: X-scores from weights
      t_new <- Xs_deflated %*% w_new
      
      # Step 3: Y-loadings from X-scores
      c_new <- t(Ys_deflated) %*% t_new / sum(t_new^2)
      
      # Step 4: Y-scores from Y-loadings
      u_new <- Ys_deflated %*% c_new / sum(c_new^2)
      
      # Convergence Test
      if (sum((w_new - w)^2) < 1e-06) break
      w <- w_new; u <- u_new
    }
    
    # Step 5: X-loadings from X-scores
    p_new <- t(Xs_deflated) %*% t_new / sum(t_new^2)
    
    # Saving components
    W[, h]  <- w_new
    Tt[, h] <- t_new
    U[, h]  <- u_new
    Cc[, h] <- c_new
    Pp[, h] <- p_new
    
    # CROSS-VALIDATION incremental 
    if (cv) {
      # Compute RSS after this component (baseline for the next iteration)
      RSS[h + 1, ] <- colSums((Ys_deflated - t_new %*% t(c_new))^2)
      
      # Cross-validation on deflated's data until h-1
      press_fold <- matrix(0, length(folds), q)
      
      for (fold_idx in seq_along(folds)) {
        idx <- folds[[fold_idx]]
        
        # Training and test on data deflated  
        Xtr <- Xs_deflated[-idx, , drop = FALSE]
        Ytr <- Ys_deflated[-idx, , drop = FALSE]
        Xte <- Xs_deflated[idx, , drop = FALSE]
        Yte <- Ys_deflated[idx, , drop = FALSE]
        
        # Algorithm PLS on training data of the fold 
        u_cv <- Ytr[, 1]; w_cv <- rep(1, p)
        for (iter in 1:100) {
          w_new_cv <- t(Xtr) %*% u_cv / sum(u_cv^2)
          w_new_cv <- w_new_cv / sqrt(sum(w_new_cv^2))
          t_new_cv <- Xtr %*% w_new_cv
          c_new_cv <- t(Ytr) %*% t_new_cv / sum(t_new_cv^2)
          u_new_cv <- Ytr %*% c_new_cv / sum(c_new_cv^2)
          if (sum((w_new_cv - w_cv)^2) < 1e-06) break
          w_cv <- w_new_cv; u_cv <- u_new_cv
        }
        
        # Predictions on tests data 
        Yhat_cv <- (Xte %*% w_new_cv) %*% t(c_new_cv)
        press_fold[fold_idx, ] <- colSums((Yte - Yhat_cv)^2)
      }
      
      PRESS[h, ] <- colSums(press_fold)
      
      # *** Q² incremental. Dynamic baseline ***
      Q2[h, ] <- 1 - (PRESS[h, ] / RSS[h, ])
    }
    
    # Deflating fot the next component 
    Xs_deflated <- Xs_deflated - t_new %*% t(p_new)
    Ys_deflated <- Ys_deflated - t_new %*% t(c_new)
  }
  
  Ws <- W %*% solve(t(Pp) %*% W, tol = 1e-20)
  sdX <- apply(X, 2, sd);  sdX[sdX == 0] <- 1
  sdY <- apply(Y, 2, sd);  sdY[sdY == 0] <- 1
  
  Bs  <- Ws %*% t(Cc)
  Br  <- sweep(Bs, 1, sdX, "/")
  Br  <- sweep(Br, 2, sdY, "*")
  cte <- colMeans(Y) - colMeans(X) %*% Br
  Yhat<- X %*% Br + matrix(rep(cte, each = n), n, q)
  
  VIP <- .vip(W, colMeans(cor(Y, Tt)^2))
  
  # Q²CUM Multiplicative 
  if (cv) {
    Q2cum <- matrix(0, nc, q)
    Q2cum[1, ] <- Q2[1, ]
    if (nc > 1) {
      for (i in 2:nc) {
        # Multiplicative Formula: Q²cum[i] = 1 - ∏(1 - Q²[j]) per j=1:i
        Q2cum[i, ] <- 1 - apply(1 - Q2[1:i, , drop=FALSE], 2, prod)
      }
    }
  } else {
    Q2 <- Q2cum <- NULL
  }
  
  structure(list(
    x.scores = Tt, x.loads = Pp, y.scores = U, y.loads = Cc,
    raw.wgs = W,   mod.wgs = Ws,
    std.coefs = Bs, reg.coefs = rbind(Br, INTERCEPT = cte),
    y.pred = Yhat, resid = Y - Yhat,
    VIP = VIP,  Q2 = Q2, Q2cum = Q2cum),
    class = "pls_reg")
}
