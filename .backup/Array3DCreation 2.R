#' Build a 3-D numeric array (subjects × features × time)
#'
#' Converts a “long” data.frame with **subject IDs**, **time points** and
#' **feature columns** into a numeric array whose dimensions are  
#' *Subject* × *Feature* × *Time*.
#'
#' @param df A data.frame in long format.
#' @param id_col   Name of the subject-identifier column.
#' @param time_col Name of the time-point column.
#' @param feature_cols Character vector of feature column names. If `NULL`,
#'   all non-`exclude_cols` are used.
#' @param exclude_cols Columns to ignore entirely.
#' @param subjects Optional character vector specifying the ordering of the
#'   *Subject* dimension.
#' @param time_points Optional numeric / character vector giving the ordering
#'   of the *Time* dimension.
#' @param legacy_na Logical; if `TRUE`, forces `deduplicate = "first"` and
#'   applies numeric coercion rules compatible with legacy code.
#' @param deduplicate Strategy when duplicate `(id_col, time_col)` pairs are
#'   found: `"first"`, `"last"`, or `"mean"`.
#'
#' @return A 3-D numeric array with dimension names *Subject*, *Feature*,
#'   and *Time*.
#'
#' @examples
#' # df must contain columns id, time, and features
#' arr <- .create_array3d(df, id_col = "ID", time_col = "T")
#' dim(arr)
#'
#' @keywords internal
#' @rdname helpers-internal

.create_array3d <- function(df, id_col, time_col,
                           feature_cols = NULL,
                           exclude_cols = c(id_col, time_col),
                           subjects     = NULL,
                           time_points  = NULL,
                           legacy_na    = FALSE,
                           deduplicate  = c("first","last","mean")) {
  deduplicate <- match.arg(deduplicate)
  if (legacy_na) deduplicate <- "first"
  
  if (is.null(feature_cols)) feature_cols <- setdiff(names(df), exclude_cols)
  if (is.null(subjects))     subjects    <- unique(df[[id_col]])
  if (is.null(time_points))  time_points <- sort(unique(df[[time_col]]))
  
  arr <- array(NA_real_,
               dim = c(length(subjects), length(feature_cols), length(time_points)),
               dimnames = list(Subject = as.character(subjects),
                               Feature = feature_cols,
                               Time    = as.character(time_points)))
  
  for (k in seq_along(time_points)) {
    tp    <- time_points[k]
    slice <- df[df[[time_col]] == tp, c(id_col, feature_cols), drop = FALSE]
    
    if (nrow(slice) > 0 && any(duplicated(slice[[id_col]]))) {
      if (deduplicate == "first") slice <- slice[!duplicated(slice[[id_col]]), ]
      if (deduplicate == "last")  slice <- slice[rev(!duplicated(rev(slice[[id_col]]))), ]
      if (deduplicate == "mean") {
        formula_str <- paste(".", "~", id_col)
        slice <- aggregate(as.formula(formula_str), data = slice, FUN = mean, na.rm = TRUE)
      }
    }
    
    anchor <- data.frame(tmp = subjects); names(anchor)[1] <- id_col
    merged <- merge(anchor, slice, by = id_col, all.x = TRUE, sort = TRUE)
    
    mat <- as.matrix(merged[, feature_cols, drop = FALSE])
    storage.mode(mat) <- "numeric"   
    arr[, , k] <- mat
  }
  arr
}
