#The function is agnostic to the VIP object name: pass `VIP3Dmodel1`,
#' `VIP3Dmodel2`, `VIP2D`, or a mixOmics `vip()` matrix.
#'This function is reccomended to be used if we have large dataset to improve NPLSDA results
#' @param vip         Matrix/data‑frame of VIP scores (variables × columns).
#' @param thr         Percentile threshold (0–100). Default 95.
#' @param cols        Columns to use (numeric indices or names). `NULL` keeps
#'                    all columns.
#' @param aggregator  How to summarise across selected columns:
#'                    * "any"  – keep if passes **one** column (logical OR)
#'                    * "mean" – percentile filter on the **row‑wise mean**
#' @param strip_time  Logical. Remove the suffix matched by `pattern`.
#' @param pattern     Regex used when `strip_time = TRUE` (default "\|t\\d+$").
#' @param na.rm       Logical passed to `quantile()` and row operations.
#'
#' @return Character vector of selected variables (deduplicated) 
#' @export
feature_selection <- function(vip,
                       thr = 95,
                       cols = NULL,
                       aggregator = c("any","mean", "all"),
                       strip_time = TRUE,
                       pattern = "\\|t\\d+$",
                       na.rm = TRUE)
                      {
  aggregator <- match.arg(aggregator)
  
  vip_df <- as.data.frame(vip)
  
  # ---- subset columns if requested
  if (!is.null(cols)) {
    if (is.numeric(cols)) {
      bad <- cols[cols < 1 | cols > ncol(vip_df)]
      if (length(bad))
        stop(sprintf("`cols` index out of range 1..%d (provided: %s)",
                     ncol(vip_df), paste(bad, collapse = ", ")), call. = FALSE)
    } else {
      missing <- setdiff(cols, colnames(vip_df))
      if (length(missing))
        stop(sprintf("`cols` has unknown column name(s): %s. Available: %s",
                     paste(missing, collapse = ", "),
                     paste(colnames(vip_df), collapse = ", ")), call. = FALSE)
    }
    vip_df <- vip_df[, cols, drop = FALSE]
  }
  
  if (ncol(vip_df) == 0)
    stop("No columns selected in vip matrix.", call. = FALSE)
  
  # ---- percentile filter column‑by‑column -----------------------------
  bool_matrix <- sapply(vip_df, function(col) {
    thr_val <- stats::quantile(col, thr / 100, na.rm = na.rm, names = FALSE)
    col >= thr_val
  })
  
  # ---- aggregate across columns
  keep <- switch(aggregator,
                 any  = apply(bool_matrix, 1, any),
                 all  = apply(bool_matrix, 1, all),
                 mean = {
                   row_mean <- rowMeans(vip_df, na.rm = na.rm)
                   # prendo per ogni riga  il thr-esimo percentile 
                   row_mean >= stats::quantile(row_mean, thr / 100, names = FALSE, na.rm = na.rm)
                 })
  
  vars <- rownames(vip_df)[keep]
  
  # ---- optional suffix removal ---------------------------------------
  if (strip_time) vars <- sub(pattern, "", vars)
  vars <- unique(vars)
  
  vars
}

# NPLSDA136GCTOFPackagePatchArrayPack2D$VIP3Dmodel2
# featureSelection(NPLSDA136GCTOFPackagePatchArrayPack2D$VIP3Dmodel1)
# featureSelection(NPLSDA136GCTOFPackagePatchArrayPack2D$VIP3Dmodel2)
# featureSelection(NPLSDA136GCTOFPackagePatchArrayPack2D$VIP2D,cols = 1,aggregator = "any")
# featureSelection(NPLSDA136GCTOFPackagePatchArrayPack2D$VIP2D,cols = 3,aggregator = "any")
# NPLSDA136GCTOFPackagePatchArrayPack2D$VIP3Dmodel2
# selectedGCTOFbydifVIPs$VIP3Dmodel2.95p.85v
# intersect(featureSelection(NPLSDA136GCTOFPackagePatchArrayPack2D$VIP3Dmodel2),selectedGCTOFbydifVIPs$VIP3Dmodel2.95p.85v))
# setdiff(selectedGCTOFbydifVIPs$VIP3Dmodel2.95p.85v,featureSelection(NPLSDA136GCTOFPackagePatchArrayPack2D$VIP3Dmodel2))                        
# NPLSDA136GCTOFPackagePatchArrayPack2D$NPLSDAQ2cum2D
# NPLSDA136GCTOF$NPLSDAQ22D

# head(NPLSDA136GCTOFPackagePatchArrayPack2D$VIP3Dmodel2)
