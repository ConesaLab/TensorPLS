#' NPLS-DA elbow analysis: cumulative Q2 logic
#'
#' @param X A 3D array of shape n x p x k (subjects x variables x time).
#' @param Y A matrix n x q or a 3D array n x q x k.
#' @param max_ncomp Integer upper bound for components. If \code{NULL},
#'   uses \code{min(10, n-1, p)}.
#' @param reps Integer, number of CV repetitions (default 1). When \code{reps > 1},
#'   the function reports standard errors of cumulative Q2 across repetitions.
#' @param seed Integer or \code{NULL}. If provided, controls fold randomization
#'   across repetitions (\code{seed + r - 1}).
#'
#' @return An object of class \code{ncomp_elbow} containing:
#' \itemize{
#'   \item \code{table}: a \code{data.frame} with columns
#'     \code{comp}, \code{Q2_inc}, \code{Q2_cum}, \code{gain_cum}, \code{Q2_cum_se}
#'   \item \code{call}: a list with the effective analysis parameters
#' }
#' \itemize{
#'   \item \code{X} is unfolded to an \code{n x (p*k)} matrix.
#'   \item \code{Y} is unfolded to an \code{n x q prime} matrix if 3D, or coerced to matrix if already 2D.
#' }
#' For each repetition, a single \code{pls_reg(Xmat, Ymat, ncomp = Hmax, cv = TRUE)}
#' run is used to extract incremental and cumulative Q2 up to H components.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' X <- array(rnorm(50*20*3), c(50, 20, 3))
#' Y <- matrix(rnorm(50*2), 50, 2)
#' res <- ncomp_elbow_nplsda(X, Y, max_ncomp = 8, reps = 1)
#' print(res)
#' plot(res)
#' head(res$table)
#' }
#'
#' @export
ncomp_elbow_nplsda <- function(X, Y, max_ncomp = NULL, reps = 1, seed = 123) {
  stopifnot(length(dim(X)) == 3L)
  n <- dim(X)[1]; p <- dim(X)[2]; k <- dim(X)[3]
  
  # Unfold X: n x (p*k)
  Xmat <- matrix(aperm(X, c(1, 2, 3)), nrow = n)
  
  # Unfold Y as in the unified workflow
  if (length(dim(Y)) == 3L) {
    Ymat <- matrix(aperm(Y, c(1, 2, 3)), nrow = n)
  } else {
    Ymat <- as.matrix(Y)
  }
  
  Hmax <- if (is.null(max_ncomp)) min(10L, n - 1L, p) else min(as.integer(max_ncomp), 10L, n - 1L, p)
  if (Hmax < 2L) stop("Not enough samples/variables to estimate >= 2 components")
  
  inc_mat <- matrix(NA_real_, nrow = Hmax, ncol = reps)
  cum_mat <- matrix(NA_real_, nrow = Hmax, ncol = reps)
  
  for (r in seq_len(reps)) {
    if (!is.null(seed)) set.seed(seed + r - 1L)
    fit <- pls_reg(Xmat, Ymat, ncomp = Hmax, cv = TRUE)  #use pls_reg
    
    inc <- rowMeans(fit$Q2,    na.rm = TRUE)  # incremental Q2 per component
    cum <- rowMeans(fit$Q2cum, na.rm = TRUE)  # cumulative Q2 per component
    
    inc_mat[, r] <- inc
    cum_mat[, r] <- cum
  }
  
  inc_mean <- rowMeans(inc_mat, na.rm = TRUE)
  cum_mean <- rowMeans(cum_mat, na.rm = TRUE)
  
  # Standard error across repetitions for cum Q2 (NA if reps == 1)
  cum_se <- if (reps > 1L) {
    apply(cum_mat, 1L, function(x) {
      sx <- sd(x, na.rm = TRUE); nx <- sum(is.finite(x))
      if (nx > 0) sx / sqrt(nx) else NA_real_
    })
  } else {
    rep(NA_real_, Hmax)
  }
  
  gains <- c(cum_mean[1L], diff(cum_mean))
  
  out_tab <- data.frame(
    comp      = seq_len(Hmax),
    Q2_inc    = inc_mean,
    Q2_cum    = cum_mean,
    gain_cum  = gains,
    Q2_cum_se = cum_se
  )
  
  structure(list(
    table = out_tab,
    call  = list(max_ncomp = Hmax, reps = reps, seed = seed)
  ), class = "ncomp_elbow")
}

#' @export
print.ncomp_elbow <- function(x, ...) {
  cat("NPLS-DA elbow object (no recommendation)\n")
  cat("Components evaluated:", nrow(x$table), "\n")
  cat("Columns:", paste(colnames(x$table), collapse = ", "), "\n")
  invisible(x)
}

#' @export
plot.ncomp_elbow <- function(x, ...) {
  op <- par(no.readonly = TRUE); on.exit(par(op))
  # Elbow plot: cumulative Q2 vs components
  plot(x$table$comp, x$table$Q2_cum, type = "b",
       xlab = "Components", ylab = "Cumulative Q2", ...)
  # Error bars if available
  if (!all(is.na(x$table$Q2_cum_se))) {
    arrows(x$table$comp, x$table$Q2_cum - x$table$Q2_cum_se,
           x$table$comp, x$table$Q2_cum + x$table$Q2_cum_se,
           angle = 90, code = 3, length = 0.05)
  }
  invisible(NULL)
}



#load("/Users/alessandrogiordano/Desktop/AnaConesa/NPLSDAPackage_Support/Imputation/Obj/fullarrayGCTOFX.RData")
#dim(fullarrayGCTOFX)
#load("/Users/alessandrogiordano/Downloads/outcomedummyarray136 (1).RData")
#suggested = ncomp_elbow_nplsda(fullarrayGCTOFX,outcomedummyarray136,reps = 10)
#suggested$table
#plot.ncomp_elbow(suggested)

