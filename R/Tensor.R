#' Prepare omics data for downstream analysis 
#'
#' Reads raw tabular data (long or wide), applies a series of
#' preprocessing steps (optional transposition, numeric coercion, cohort
#' filtering, duplicate handling, etc.) and returns either
#' *-* a three-dimensional **tensor**
#'    (*Subject × Feature × Time*) *or*
#' *-* a two-dimensional **matrix**
#'    (*Subject × Feature*), depending on the presence of a time column or
#'   the `tensor` flag.
#'
#' @section Workflow in brief:
#' \enumerate{
#'   \item \strong{Input loading}\\
#'         `data` and, if applicable, `cohort` may be in-memory objects or
#'         file paths (`csv`, `tsv`, `rds`, `feather`, `parquet`).
#'   \item \strong{Transpose (optional)} via \code{.transpose_if_needed()}.
#'   \item \strong{Numeric coercion} (legacy vs.\ modern rules).
#'   \item \strong{Cohort filter} applied on an auxiliary table.
#'   \item \strong{Minimum time-points} filter.
#'   \item \strong{Tensor / matrix creation} through
#'         \code{.build_tensor()}.
#' }
#'
#' @param data Data.frame, matrix/array, or a readable file path.
#' @param id_col Character; name of the subject identifier column.
#' @param time_col Character; name of the time-point column. If `NULL`
#'   or not present in `data`, a 2-D matrix is produced.
#' @param transpose One of `"auto"`, `"never"`, `"always"`; passed to
#'   \code{.transpose_if_needed()}.
#' @param header_in_row Logical; used only when transposing.
#' @param cohort Optional data.frame or file path containing metadata.
#' @param cohort_id_col Column in `cohort` matching `id_col`.
#' @param cohort_filter Character scalar; an R expression evaluated within
#'   `cohort` to subset rows.
#' @param min_timepoints Integer; keep only subjects observed at least this
#'   many times.
#' @param tensor Logical; if `FALSE`, return the (filtered) data.frame
#'   rather than tensor/matrix.
#' @param exclude_cols Columns to drop entirely when building features.
#' @param feature_cols Character vector; explicit feature columns to keep.
#'   Overrides `exclude_cols`.
#' @param subjects Optional character vector defining the order (and subset)
#'   of subjects in the tensor/matrix.
#' @param time_points Optional vector defining the order (and subset) of
#'   time points in the tensor.
#' @param numeric_coercion Logical; if `TRUE` attempt to coerce all columns
#'   to numeric (ignored when `legacy_na = TRUE`).
#' @param legacy_na Logical; reproduces historical NA-handling rules and
#'   forces \code{deduplicate = "first"}.
#' @param deduplicate Strategy used when duplicate subject/time pairs are
#'   encountered: `"first"`, `"last"`, or `"mean"`.
#'
#' @return
#' \itemize{
#'   \item \strong{3-D numeric array} (Subj × Feat × Time) if
#'         \code{tensor = TRUE} and a valid `time_col` is given;
#'   \item \strong{2-D numeric matrix} (Subj × Feat) if
#'         \code{tensor = TRUE} but no time column is available;
#'   \item \strong{data.frame} after filtering/coercion if
#'         \code{tensor = FALSE}.
#' }
#'
#'
#' @examples
#' ## Minimal runnable example using package-internal demo files
#' raw_path    <- system.file("extdata", "GCTOF_Data_Processed.csv",
#'                            package = "NPLS-Multiomi")
#' cohort_path <- system.file("extdata", "CohortData.csv",
#'                            package = "NPLS-Multiomi")
#'
#' arr <- prepare_omics(
#'   data          = raw_path,
#'   id_col        = "Individual.Id",
#'   time_col      = "Time.to.IA",
#'   transpose     = "always",
#'   legacy_na     = TRUE,
#'   cohort        = cohort_path,
#'   cohort_id_col = "Group.Id",
#'   cohort_filter = "Model.or.Validation=='Model'"
#' )
#' dim(arr)
#'
#' @export
#' @importFrom yaml read_yaml
#' @importFrom arrow read_feather read_parquet
#' 
#' 
#' 
prepare_omics <- function(
    data,
    id_col,
    time_col           = NULL,
    transpose          = "auto",
    header_in_row      = TRUE,
    cohort             = NULL,
    cohort_id_col      = id_col,
    cohort_filter      = NULL,
    min_timepoints     = NULL,
    tensor             = TRUE,
    exclude_cols       = NULL,
    feature_cols       = NULL,
    subjects           = NULL,
    time_points        = NULL,
    numeric_coercion   = TRUE,
    legacy_na          = FALSE,
    deduplicate        = c("first","last","mean")) {
  

  deduplicate <- match.arg(deduplicate)
  
  # About I/O & transpose
  obj <- .read_or_pass(data)
  if (is.matrix(obj) || (is.array(obj) && length(dim(obj)) == 3L)) return(obj)
  if (!is.data.frame(obj)) stop("Unsupported input type: ", class(obj))
  
  obj <- .transpose_if_needed(obj, transpose, header_in_row)
  
  # About (2) Numeric coercion
  if (legacy_na) {
    obj[[id_col]] <- suppressWarnings(as.numeric(as.character(obj[[id_col]])))
    if (!is.null(time_col) && time_col %in% names(obj))
      obj[[time_col]] <- suppressWarnings(as.numeric(as.character(obj[[time_col]])))
  } else if (numeric_coercion) {
    obj[] <- lapply(obj, function(x) suppressWarnings(as.numeric(as.character(x))))
  }
  
  # About (3) Cohort filter
  if (!is.null(cohort)) {
    cdf <- .read_or_pass(cohort)
    if (!is.data.frame(cdf)) stop("Cohort must be data.frame or file.")
    if (cohort_id_col %in% names(cdf))
      names(cdf)[names(cdf) == cohort_id_col] <- id_col
    if (!is.null(cohort_filter))
      cdf <- cdf[with(cdf, eval(parse(text = cohort_filter))), ]
    obj <- obj[obj[[id_col]] %in% cdf[[id_col]], ]
  }
  
  # About (4) min_timepoints 
  if (!is.null(min_timepoints) && !is.null(time_col) && time_col %in% names(obj)) {
    keep <- names(which(table(obj[[id_col]]) >= min_timepoints))
    obj  <- obj[obj[[id_col]] %in% keep, ]
  }
  
  # (5) Tensor / Matrix
  if (!tensor) return(obj)
  
  .build_tensor(obj,
                id_col       = id_col,
                time_col     = time_col,
                exclude_cols = if (is.null(exclude_cols)) c(id_col, time_col) else exclude_cols,
                feature_cols = feature_cols,
                subjects     = subjects,
                time_points  = time_points,
                legacy_na    = legacy_na,
                deduplicate  = deduplicate)
}

