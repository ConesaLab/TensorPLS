#' Imputemethod  –  Tucker-3 iterative imputation (Imputation of missing values via iterative Tucker-3)
#'
#' @param X        3-way array (subjects × variables × time)
#' @param fac      integer length-3: n° components for each mode (default c(2,2,2))
#' @param conver   defauld threshold convergence ΔSSTc (default 1e-7)
#' @param max.iter max number of iteration (default 1000)
#' @param seed     integer
#' @param verbose  default(TRUE)
#'
#' @return The array \code{X} with \code{NA} replaced
#' @export
Imputemethod <- function (X, fac = c(2, 2, 2),
                                 conver   = 1e-7,
                                 max.iter = 1000,
                                 seed     = NULL,
                                 verbose  = TRUE)
{
  ##1st check 
  stopifnot(is.array(X),
            length(dim(X)) == 3L,
            is.numeric(fac),
            length(fac) == 3L)

  if (!is.null(seed)) set.seed(seed)

  if (!anyNA(X)) {
    if (verbose) message("No NA values present – imputation skipped.")
    return(X)
  }

  ## Initialization of NA using rnorm
  if (verbose) message("Imputing missing values with iterative Tucker-3 …")
  d        <- dim(X)
  NA.pos   <- which(is.na(X))
  mu	   <- mean(X, na.rm = TRUE)
  sigma    <- sd(X,   na.rm = TRUE)
  if (is.na(sigma) || sigma == 0) sigma <- 1           # edge-case
  X[NA.pos] <- rnorm(length(NA.pos), mu, sigma)

  ##Factors clipping to avoid numbers of factors > size array's size
  fac <- pmin(fac, d)

  ## Tucker 3 Iterations
  SSTc.old <- 0
  for (it in seq_len(max.iter)) {

    Xe <- tucker3mod2(X,
                             COMP     = fac,
                             conver   = 1e-7,
                             max.iter = 500,
                             seed     = seed,
                             verbose  = FALSE)

    X[NA.pos] <- Xe$Xhat[NA.pos]

    SSTc      <- sum(Xe$G^2)
    delta     <- abs(SSTc - SSTc.old) / SSTc
    if (delta < conver) break
    SSTc.old  <- SSTc
  }

  if (it == max.iter && verbose)
    message("Reached max.iter (", max.iter, ") – delta = ", formatC(delta))

  X
}
 
#' Tucker-3 decomposition of a 3-way array.
#'
#' @param X       3-way array (subjects × variables × time-points)
#' @param COMP    integer length-3: number of components (P, Q, R)
#' @param conver  numeric, relative threshold for convergence on SSF
#' @param max.iter integer, maximum ALS iterations
#' @param seed    optional integer for reproducible random start
#' @param verbose logical; print progress if TRUE
#'
#' @return list with Factors, core array, Xhat, SST, SSF, SSE, expl.var,
#'         GCV, edf and tdf – identical structure to the original author.
#' @export
tucker3mod2 <- function (X, COMP,
                                conver   = 1e-7,
                                max.iter = 10000,
                                seed     = NULL,
                                verbose  = FALSE)
{
  ## -- sanity-check
  stopifnot(is.array(X), length(dim(X)) == 3L,
            is.numeric(COMP), length(COMP) == 3L)

  if (any(COMP > dim(X)))
    stop("COMP components exceed array dimensions: ",
         paste(COMP, collapse = ","), " vs ",
         paste(dim(X), collapse = ","))

  if (!is.null(seed)) set.seed(seed)

  ## -- imputation fallback
  if (anyNA(X)) {
    if (verbose) message("Missing values detected – imputing internally …")
    X <- Imputemethod(X,
                             fac  = pmin(COMP, dim(X)),
                             seed = seed,
                             verbose = FALSE)
  }

  d  <- dim(X)                                   # c(N, P, T)
  XA <- matrix(X,           d[1], d[2] * d[3])   # unfold mode-1
  XB <- matrix(aperm(X, c(2,1,3)), d[2], d[1] * d[3])   # mode-2
  XC <- matrix(aperm(X, c(3,1,2)), d[3], d[1] * d[2])   # mode-3

  ## -- random orthonormal starts (SVD on noise) -----------------------------
  A <- svd(matrix(rnorm(d[1] * COMP[1]), d[1], COMP[1]), nu = COMP[1])$u
  B <- svd(matrix(rnorm(d[2] * COMP[2]), d[2], COMP[2]), nu = COMP[2])$u
  C <- svd(matrix(rnorm(d[3] * COMP[3]), d[3], COMP[3]), nu = COMP[3])$u

  SST <- SSF <- SSF.old <- 0
  if (verbose)
    message("iter   SST            SSF            SSE         expl.var")

  ## -- ALS loop
  for (it in seq_len(max.iter)) {

    ## update factors
    A1 <- svd(XA %*% kronecker(C,  B), nu = COMP[1])$u
    B1 <- svd(XB %*% kronecker(C,  A), nu = COMP[2])$u
    C1 <- svd(XC %*% kronecker(B,  A), nu = COMP[3])$u

    ## core & reconstruction
    Gu      <- t(A1) %*% XA %*% kronecker(C1, B1)
    Xhat.u  <- A1 %*% Gu %*% kronecker(t(C1), t(B1))

    ## fit statistics
    SST <- sum(X^2)
    SSF <- sum(Xhat.u^2)
    SSE <- SST - sum(rowSums(Gu^2))
    expl <- SSF / sum(XA^2)

    if (verbose && it %% 50 == 1)
      message(sprintf("%4d  %.3e  %.3e  %.3e  %.5f",
                      it, SST, SSF, SSE, expl))

    ## ---- convergence test
    if (abs(SSF - SSF.old) / SST < conver) break
    SSF.old <- SSF                                    

    if (sum((A1 - A)^2) < conver &&
        sum((B1 - B)^2) < conver &&
        sum((C1 - C)^2) < conver) break

    ## prepare next iteration
    A <- A1; B <- B1; C <- C1
  }

  if (it == max.iter && verbose)
    message("Maximum iterations (", max.iter, ") reached – last ΔSSF/SST = ",
            signif(abs(SSF - SSF.old) / SST, 3))

  ## -- wrap up-
  G    <- array(as.vector(Gu), COMP)
  Xhat <- array(as.vector(Xhat.u), dim = d, dimnames = dimnames(X))

  edf <- c(d[1]*COMP[1] - COMP[1]*(COMP[1]+1)/2,
           d[2]*COMP[2] - COMP[2]*(COMP[2]+1)/2,
           d[3]*COMP[3] - COMP[3]*(COMP[3]+1)/2,
           prod(COMP))
  tdf     <- sum(d*COMP) + sum(COMP) - sum(COMP^2)
  pxdim   <- prod(d)
  GCV     <- (SSE/pxdim) / (1 - sum(edf)/pxdim)^2

  ## dim-names for factors
  rownames(A1) <- dimnames(X)[[1]]
  rownames(B1) <- dimnames(X)[[2]]
  rownames(C1) <- dimnames(X)[[3]]

  list(FactorsX = list(Mode1 = A1,
                       Mode2 = B1,
                       Mode3 = C1),
       G         = G,
       Xhat	 = Xhat,
       SST	 = SST,
       SSF	 = SSF,
       SSE	 = SSE,
       expl.var  = expl,
       GCV	 = GCV,
       edf	 = edf,
       tdf	 = tdf)
}
