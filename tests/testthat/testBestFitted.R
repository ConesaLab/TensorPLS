# Config/testthat/edition: 3

library(testthat)

# Small Array 
make_tiny_array <- function() {
  set.seed(123)
  array(rnorm(4 * 3 * 2), dim = c(4, 3, 2))
}

# Mock: centering -> identity. No centering
fake_center_mode <- function(X, mode) X

# Mock: imputation -> identity
fake_impute <- function(Xc) Xc

# Mock: Tucker3.
# - For a grid P,Q,R in {1,2}, we assign defalut EXP.VAR
# - Checked failed for (2,1,1) to test handling errors.
fake_tucker <- function(A3D, COMP, conver = 1e-6, max.iter = 500, verbose = FALSE) {
  P <- as.integer(COMP[1]); Q <- as.integer(COMP[2]); R <- as.integer(COMP[3])
  key <- paste(P, Q, R, sep = "-")
  
  # Failing for a combination 
  if (key == "2-1-1") stop("synthetic fit failure")
  
  expl_map <- list(
    "1-1-1" = 0.40,
    "1-1-2" = 0.50,
    "1-2-1" = 0.60,
    "2-2-1" = 0.70,
    "2-1-2" = 0.65,
    "1-2-2" = 0.68,
    "2-2-2" = 0.75
  )
  expl <- expl_map[[key]]
  if (is.null(expl)) expl <- 0.10
  
  list(expl.var = expl, SSE = 1 - expl, tdf = P + Q + R, SSF = expl)
}


# Test
test_that("bestfittedmodel returns expected structure and columns", {
  X <- make_tiny_array()
  
  testthat::local_mocked_bindings(
    center_mode           = fake_center_mode,
    Imputemethod   = fake_impute,
    tucker3mod2    = fake_tucker
  )
  
  # Small grid : P,Q,R in {1,2} -> 2^3 = 8 models (1 models go to fails but expect to see the row with na)
  out <- bestfittedmodel(X, centering = 2, min_comp = 1, max_comp = 2, verbose = FALSE)
  
  # Type and classes
  expect_type(out, "list")
  expect_true("bfm" %in% class(out))
  
  # Names provided by the output
  expect_setequal(names(out), c("summary", "top_by_fit", "note"))
  
  # Columns expected 
  expect_equal(nrow(out$summary), 8)
  expect_true(all(c("P", "Q", "R", "EXPL.VAR", "SSE", "DF", "SSF", "SumComp", "Fitper") %in%
                    colnames(out$summary)))
  
  # Fitper = EXPL.VAR * 100 
  non_na <- subset(out$summary, !is.na(EXPL.VAR))
  expect_equal(non_na$Fitper, non_na$EXPL.VAR * 100)
  
  # Row that has failed fit should have NA 
  expect_equal(sum(is.na(out$summary$Fitper)), 1)
  
  #Note
  expect_true(grepl("^Returned full grid summary", out$note))
})


# Test: ordinamento "parsimonia-first" poi Fit%

test_that("rows are ordered by SumComp ascending, then Fit% descending (among non-NA)", {
  X <- make_tiny_array()
  
  testthat::local_mocked_bindings(
    center_mode           = fake_center_mode,
    Imputemethod   = fake_impute,
    tucker3mod2    = fake_tucker
  )
  
  out <- bestfittedmodel(X, centering = 2, min_comp = 1, max_comp = 2, verbose = FALSE)
  
  # Consideriamo solo le righe con Fitper non-NA per verificare l'ordinamento
  tab <- subset(out$summary, !is.na(Fitper))
  
  # 1) SumComp no descending
  expect_true(all(diff(tab$SumComp) >= 0))
  
  # 2) In SumComp==5 we expect this order of Fit%: 70 > 68 > 65
  g5 <- subset(tab, SumComp == 5)
  expect_equal(g5$Fitper, sort(g5$Fitper, decreasing = TRUE))
})

# Test: top_by_fit restituisce i modelli con Fit% maggiore (ignorando NA)
test_that("top_by_fit returns highest-Fit% rows with proper rownames", {
  X <- make_tiny_array()
  
  testthat::local_mocked_bindings(
    center_mode           = fake_center_mode,
    Imputemethod   = fake_impute,
    tucker3mod2    = fake_tucker
  )
  
  out <- bestfittedmodel(X, centering = 2, min_comp = 1, max_comp = 2,
                         top_fit_n = 3, verbose = FALSE)
  
  # Expecting 3x4 table (P,Q,R,Fitper)
  expect_s3_class(out$top_by_fit, "data.frame")
  expect_equal(nrow(out$top_by_fit), 3)
  expect_setequal(colnames(out$top_by_fit), c("P","Q","R","Fitper"))
  
  # Ordered by descending 
  expect_equal(out$top_by_fit$Fitper, sort(out$top_by_fit$Fitper, decreasing = TRUE))
  
  # 3 best 75, 70, 68
  expect_equal(round(out$top_by_fit$Fitper, 2), c(75, 70, 68))
  
  # Rownames "TopFit1", "TopFit2", ...
  expect_true(all(startsWith(rownames(out$top_by_fit), "TopFit")))
})

# Test: verbose stampa un messaggio con numero di modelli
test_that("verbose prints a brief summary message", {
  X <- make_tiny_array()
  
  testthat::local_mocked_bindings(
    center_mode           = fake_center_mode,
    Imputemethod   = fake_impute,
    tucker3mod2    = fake_tucker
  )
  
  # 8 modelli in  the grid (one failed but the row is stil present)
  expect_message(
    bestfittedmodel(X, centering = 2, min_comp = 1, max_comp = 2,
                    top_fit_n = 2, verbose = TRUE),
    regexp = "Full grid evaluated: 8 models\\."
  )
})

# Test: input non valido per 'centering'
test_that("invalid centering throws a clear error", {
  X <- make_tiny_array()
  
  expect_error(
    bestfittedmodel(X, centering = 9, min_comp = 1, max_comp = 2, verbose = FALSE),
    regexp = "centering must be 0,1,2,3"
  )
})

# Test: top_fit_n NULL -> top_by_fit è NULL
test_that("top_fit_n = NULL yields top_by_fit = NULL", {
  X <- make_tiny_array()
  
  testthat::local_mocked_bindings(
    center_mode           = fake_center_mode,
    Imputemethod   = fake_impute,
    tucker3mod2    = fake_tucker
  )
  
  out <- bestfittedmodel(X, centering = 0, min_comp = 1, max_comp = 2,
                         top_fit_n = NULL, verbose = FALSE)
  expect_null(out$top_by_fit)
})
