#' Build a 3-D numeric array (subjects × features × time) -- robust version
#' 
#' @keywords internal
.create_array3d <- function(df, id_col, time_col,
                            feature_cols = NULL,
                            exclude_cols = c(id_col, time_col),
                            subjects     = NULL,
                            time_points  = NULL,
                            legacy_na    = FALSE,
                            deduplicate  = c("first","last","mean")) {
  deduplicate <- match.arg(deduplicate)
  
  # -- helpers 
  .num <- function(x) {
    # conversione numerica "tollerante": gestisce virgole decimali
    suppressWarnings(as.numeric(gsub(",", ".", as.character(x), fixed = FALSE)))
  }
  .is_date_time <- function(x) inherits(x, c("Date","POSIXct","POSIXt"))
  
  # -- validazioni minime 
  if (!is.data.frame(df)) stop("'.create_array3d' richiede un data.frame.")
  if (!all(c(id_col, time_col) %in% names(df))) {
    stop("Colonne chiave mancanti: attese '", id_col, "' e '", time_col, "'.")
  }
  
  # -- normalizza ID come character
  df[[id_col]] <- as.character(df[[id_col]])
  
  # -- determina feature
  if (is.null(feature_cols)) {
    feature_cols <- setdiff(names(df), exclude_cols)
  }
  if (length(feature_cols) == 0L) {
    stop("Nessuna feature disponibile dopo 'exclude_cols' / 'feature_cols'.")
  }
  if (!all(feature_cols %in% names(df))) {
    miss <- setdiff(feature_cols, names(df))
    stop("Le seguenti feature non esistono in 'df': ", paste(miss, collapse = ", "))
  }
  
  # -- normalizza e definisci time_points 
  
  tc <- df[[time_col]]
  
  # Teniamo traccia della "classe" target per il confronto
  time_is_date <- .is_date_time(tc)
  
  if (time_is_date) {
    # Lascia il tipo come Date/POSIX*; converti eventualmente time_points
    if (is.null(time_points)) {
      time_points <- sort(unique(tc))
    } 
    else {
      if (inherits(tc, "Date")) {
        time_points <- as.Date(time_points)
      } 
      else {
        # se tc è POSIXct, usa stessa timezone se disponibile
        tz <- attr(tc, "tzone"); if (is.null(tz)) tz <- ""
        time_points <- as.POSIXct(time_points, tz = tz)
      }
    }
  } 
  else {
    # non è una data: valuta se è "numeric-like" o stringa semplice 
    tc_as_char <- as.character(tc)
    # +4 -4 " +4 " " -4 "
    is_num_like <- grepl("^\\s*[-+]?\\d*(?:[\\.,]\\d+)?\\s*$", tc_as_char) | is.na(tc)
    if (all(is_num_like, na.rm = TRUE)) {
      df[[time_col]] <- .num(tc_as_char)
      if (is.null(time_points)) {
        time_points <- sort(unique(df[[time_col]]))
      } else {
        time_points <- .num(time_points)
      }
      time_is_date <- FALSE
    } 
    else {
      # NO num. NO Date --> is String
      df[[time_col]] <- tc_as_char
      if (is.null(time_points)) {
        # ordine di apparizione
        time_points <- unique(df[[time_col]])
      } else {
        time_points <- as.character(time_points)
      }
    }
  }
  
  # -- definisci subjects (ordine/insieme) 
  if (is.null(subjects)) {
    subjects <- unique(df[[id_col]])
  }
  subjects <- as.character(subjects)
  
  # -- prepara array destinazione 
  arr <- array(NA_real_,
               dim = c(length(subjects), length(feature_cols), length(time_points)),
               dimnames = list(Subject = subjects,
                               Feature = feature_cols,
                               Time    = as.character(time_points)))
  # Nota: dimnames$Time viene reso character per coerenza con il resto del pacchetto
  
  # -- ciclo sui time point
  for (k in seq_along(time_points)) {
    tp <- time_points[k]
    
    # confronto coerente col tipo selezionato
    if (time_is_date) {
      # Date/POSIXct: confronto diretto
      slice <- df[df[[time_col]] == tp, c(id_col, feature_cols), drop = FALSE]
    } 
    else if (is.numeric(df[[time_col]])) {
      slice <- df[df[[time_col]] == as.numeric(tp), c(id_col, feature_cols), drop = FALSE]
    } 
    else {
      slice <- df[df[[time_col]] == as.character(tp), c(id_col, feature_cols), drop = FALSE]
    }
    
    # assicurati che la chiave sia character
    if (nrow(slice) > 0L) {
      slice[[id_col]] <- as.character(slice[[id_col]])
    }
    
    #  gestisci duplicati per subject
    if (nrow(slice) > 0L && any(duplicated(slice[[id_col]]))) {
      if (deduplicate == "first") {
        slice <- slice[!duplicated(slice[[id_col]]), , drop = FALSE]
      } 
      else if (deduplicate == "last") {
        slice <- slice[!duplicated(slice[[id_col]], fromLast = TRUE), , drop = FALSE]
      }
      else if (deduplicate == "mean") {
        # aggrega solo le feature, raggruppando per id_col
        agg <- aggregate(slice[feature_cols],
                         by = list(.id = slice[[id_col]]),
                         FUN = mean, na.rm = TRUE)
        names(agg)[1] <- id_col
        slice <- agg
      }
    }
    
    # -- costruisci anchor e merge preservando l'ordine dei subjects ---------
    anchor <- data.frame(tmp = subjects, stringsAsFactors = FALSE)
    names(anchor)[1] <- id_col
    anchor[[id_col]] <- as.character(anchor[[id_col]])
    
    merged <- merge(anchor, slice, by = id_col, all.x = TRUE, sort = FALSE)
    
    # riordino "blindato" per seguire esattamente 'subjects'
    ord <- match(anchor[[id_col]], merged[[id_col]])
    # se qualche subject non matcha (NA), lasciamo NA nelle feature
    if (anyNA(ord)) {
      # ricostruiamo skeleton vuoto per i soggetti mancanti
      missing_rows <- is.na(ord)
      fill <- matrix(NA_real_, nrow = sum(missing_rows), ncol = length(feature_cols))
      colnames(fill) <- feature_cols
      # righe valide
      merged_valid <- merged[!is.na(ord), , drop = FALSE]
      # costruisci merged riordinato
      merged <- rbind(
        merged_valid[match(anchor[[id_col]][!missing_rows],
                           merged_valid[[id_col]]), , drop = FALSE],
        cbind(!!id_col := anchor[[id_col]][missing_rows], as.data.frame(fill))
      )
      names(merged)[1] <- id_col
      # riordina di nuovo secondo anchor completo
      ord2 <- match(anchor[[id_col]], merged[[id_col]])
      merged <- merged[ord2, , drop = FALSE]
    } 
    else {
      merged <- merged[ord, , drop = FALSE]
    }
    
    # -- conversione esplicita delle feature a numerico 
    if (nrow(merged) > 0L) {
      for (cc in feature_cols) {
        merged[[cc]] <- .num(merged[[cc]])
      }
    }
    
    # -- estrai matrice e assegna allo slice 
    mat <- as.matrix(merged[, feature_cols, drop = FALSE])
    storage.mode(mat) <- "double"
    arr[, , k] <- mat
  }
  
  arr
}
