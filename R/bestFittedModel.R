#' Grid search of Tucker-3 models over (P, Q, R)
#'
#' Fits a Tucker-3 decomposition for all component triplets \code{(P, Q, R)}
#' within the given integer ranges and returns a complete summary table.
#'
#' 
#' The function enumerates all triplets \code{P, Q, R} in
#' \code{min_comp:max_comp} and, for each candidate, runs
#' \code{tucker3mod2()} on a preliminarily centered and imputed array.
#' Results are collated into a single data frame including explained variance
#' and error metrics. Rows are ordered by increasing model complexity
#' (\code{SumComp = P + Q + R}) and then by decreasing explained variance.
#'
#' A simple informative quantity \code{CriVal} is also returned
#' (\emph{not} a decision threshold): \code{norm(mat, "O") / (Smax - 3)},
#' where \code{mat} is a mode-1 unfolded version of the imputed tensor and
#' \code{Smax = max(P+Q+R)}.
#'
#' Optionally, the \code{top_by_fit} element contains the top-\code{N} rows
#' with highest explained variance across the whole grid.
#'
#' @param X A 3D numeric array (subjects × variables × time or any 3-way layout).
#' @param centering Integer in \code{0,1,2,3}. Which mode to center:
#'   \code{0} = no centering; \code{1,2,3} = center along the corresponding mode.
#'   Delegated to \code{center_mode()}.
#' @param min_comp Integer (>= 1). Lower bound for \code{P, Q, R}.
#' @param max_comp Integer (>= min_comp). Upper bound for \code{P, Q, R}.
#' @param top_fit_n Optional integer. If > 0, also return the top-N triplets
#'   by explained variance as \code{$top_by_fit}.
#' @param verbose Logical. If \code{TRUE}, prints a brief textual summary.
#'
#' @return A list with the following elements:
#' \item{summary}{A \code{data.frame} with one row per triplet and columns:
#'   \code{P, Q, R, EXPL.VAR, SSE, DF, SSF, SumComp, Fitper}.}
#' \item{top_by_fit}{\code{data.frame} of top-N rows by \code{Fitper}, or \code{NULL}.}
#' \item{note}{Character note describing what was returned.}
#'
#' @section Columns in \code{$summary}:
#' \describe{
#'   \item{P, Q, R}{Component counts along the three modes.}
#'   \item{EXPL.VAR}{Explained variance (0–1).}
#'   \item{SSE}{Sum of squared errors.}
#'   \item{DF}{Total degrees of freedom used by the model.}
#'   \item{SSF}{Sum of squares fitted.}
#'   \item{SumComp}{\code{P + Q + R}.}
#'   \item{Fitper}{Explained variance in percent: \code{EXPL.VAR * 100}.}
#' }
#'
#' @seealso \code{\link{plot_bfm_pareto}}, \code{\link{plot_bfm_heatmap}}
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- array(rnorm(20*15*6), c(20,15,6))
#' out <- bestfittedmodel(X, centering = 2, min_comp = 2, max_comp = 4,
#'                               top_fit_n = 5, verbose = TRUE)
#' head(out$summary)
#' }
#' @export
bestfittedmodel <- function(X, centering = 2,
                                   min_comp  = 2,
                                   max_comp  = 4,
                                   top_fit_n = NULL,
                                   verbose   = TRUE) {

  ## 0) Centering -------------------------------------------------------------
  Xc <- switch(as.character(centering),
               "0" = X,
               "1" = center_mode(X, 1),
               "2" = center_mode(X, 2),
               "3" = center_mode(X, 3),
               stop("centering must be 0,1,2,3"))

  ## 1) Preliminary imputation -----------------------------------------------
  A3D  <- Imputemethod(Xc)
  xdim <- dim(A3D)

  ## 2) Enumerate all triplets ------------------------------------------------
  grid <- expand.grid(P = min_comp:max_comp,
                      Q = min_comp:max_comp,
                      R = min_comp:max_comp)

  res <- lapply(seq_len(nrow(grid)), function(i) {
    pqr <- unlist(grid[i, ])
    mod <- try(
      tucker3mod2(A3D, COMP = pqr,
                         conver = 1e-6, max.iter = 500, verbose = FALSE),
      silent = TRUE)
    if (inherits(mod, "try-error")) {
      data.frame(P = pqr[1], Q = pqr[2], R = pqr[3],
                 EXPL.VAR = NA_real_, SSE = NA_real_,
                 DF = NA_real_, SSF = NA_real_)
    } else {
      data.frame(P = pqr[1], Q = pqr[2], R = pqr[3],
                 EXPL.VAR = mod$expl.var,
                 SSE      = mod$SSE,
                 DF       = mod$tdf,
                 SSF      = mod$SSF)
    }
  })

  tab <- do.call(rbind, res)
  tab$SumComp <- tab$P + tab$Q + tab$R
  tab$Fitper  <- tab$EXPL.VAR * 100

  ## 3) Order: parsimony first, then explained variance ----------------------
  o <- order(tab$SumComp, -tab$Fitper, na.last = TRUE)
  tab <- tab[o, ]

  ## 4) Informative critical value -------------------------------------------
  mat    <- matrix(A3D, xdim[1], xdim[2] * xdim[3])
  Smax   <- max(tab$SumComp, na.rm = TRUE)

  ## 5) Optional: top-by-fit table -------------------------------------------
  non_na <- subset(tab, !is.na(Fitper))
  if (!is.null(top_fit_n) && top_fit_n > 0 && nrow(non_na) > 0) {
    top_by_fit <- head(non_na[order(-non_na$Fitper),
                              c("P","Q","R","Fitper")], top_fit_n)
    rownames(top_by_fit) <- paste0("TopFit", seq_len(nrow(top_by_fit)))
  } else {
    top_by_fit <- NULL
  }

  note <- "Returned full grid summary. Choose about a balance between parsimony and Fit Percentage. See plots."

  if (verbose) {
    message("Full grid evaluated: ", nrow(tab), " models.")
    if (!is.null(top_by_fit)) {
      message("Top-", top_fit_n, " by Fit% (descriptive):")
      print(top_by_fit)
    }
  }

  out <- list(
    summary    = tab,
    top_by_fit = top_by_fit,
    note       = note
  )
  class(out) <- c("bfm", class(out))
  out
}



#Function for centering
# =====================================================================
#  mode = 0  → No centering
#  mode = 1  → mean for 1st mode       (rows)
#  mode = 2  → mean for 2nd mode       (columns)
#  mode = 3  → mean for 3rd mode      

center_mode <- function(X, mode = 2) {
  stopifnot(length(dim(X)) == 3L, mode %in% 0:3)
  if (mode == 0) return(X)

  d  <- dim(X)
  dn <- dimnames(X)
  perm <- switch(as.character(mode),
                 "1" = c(1, 2, 3),
                 "2" = c(2, 1, 3),
                 "3" = c(3, 1, 2))
  mat <- matrix(aperm(X, perm), d[perm[1]], prod(d[-perm[1]]))
  storage.mode(mat) <- "double"  
  
  mat <- mat - matrix(rowMeans(mat, na.rm = TRUE), nrow(mat), ncol(mat))

  Xc  <- aperm(array(mat, d[perm]), order(perm))         # back-permute
  dimnames(Xc) <- dn                                     
  Xc
}

