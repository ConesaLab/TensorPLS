#' Plot Top-N VIP2D with Case/Control strip and permutation p-values 
#'
#
#'
plot_vip2d_with_groups_nogaps <- function(
    obj, X, groups,
    comp = 1, top_n = 15, threshold = 1, sep = "_",
    mode = c("winner","effect"),
    case_level = NULL,
    case_col = "#D55E00", control_col = "#009E73",
    loser_alpha = 0.25,
    # --- permutation p-values options ---
    Y = NULL, ncomp = NULL, permute_R = 0, alpha_perm = 0.05,
    fit_fun = NULL, seed = NULL, workers = NULL,        
    # CRAN-safe options
    cap_blas_threads = FALSE,                             
    parallel = c("auto","off"),                           
    future_packages = NULL, 
    ...
) {
  mode <- match.arg(mode)
  parallel <- match.arg(parallel)
  if (!is.null(workers)) {                               
    message("Note: 'workers' is ignored; set future::plan() outside if you want parallelism.")
  }
  
  #  small helpers
  .split_feature_time <- function(keys, sep = "_") {
    feat <- sub(paste0("^(.*)", sep, "(.*)$"), "\\1", keys)
    tim  <- sub(paste0("^(.*)", sep, "(.*)$"), "\\2", keys)
    list(feature = feat, time = tim)
  }
  .coerce_groups <- function(groups, sample_ids, case_level = NULL) {
    if (is.data.frame(groups)) {
      ids  <- as.character(groups[[1]]); vals <- groups[[2]]
      idx <- match(as.character(sample_ids), ids)
      if (anyNA(idx)) stop("Some sample IDs in X not found in groups data.frame.")
      g <- vals[idx]
    } else if (!is.null(names(groups))) {
      idx <- match(as.character(sample_ids), names(groups))
      if (anyNA(idx)) stop("Some sample IDs in X not found in names(groups).")
      g <- groups[idx]
    } else {
      if (length(groups) != length(sample_ids))
        stop("groups length (", length(groups), ") != nrow(X) (", length(sample_ids), ").")
      g <- groups
    }
    g <- as.factor(g)
    if (nlevels(g) != 2)
      stop("`groups` must have exactly two levels; got: ", paste(levels(g), collapse = ", "))
    levs <- levels(g)
    if (!is.null(case_level)) {
      if (!case_level %in% levs) stop("case_level '", case_level, "' not in: ", paste(levs, collapse=", "))
      g <- factor(g, levels = c(setdiff(levs, case_level)[1], case_level))  # Control, Case
    } else if ("Class1" %in% levs) {
      g <- factor(g, levels = c(setdiff(levs, "Class1")[1], "Class1"))
    }
    g
  }
  
  .fit_fun_default <- function(Xb, Yb, ncomp, ...) {
    objb <- nplsda_vips(Xb, Yb, ncomp = ncomp, ...)
    objb$VIP2D
  }
  
  # return only intersected Top-N hits to save memory
  .permute_top_hits_once <- function(X, Y, comp, top_n, ncomp, sep, fit_fun, rn_target, perm_idx, ...) {
    Yp <- if (length(dim(Y))==3) Y[perm_idx, , , drop = FALSE] else Y[perm_idx, , drop = FALSE]
    VIP2D <- fit_fun(X, Yp, ncomp = ncomp, ...)
    rn <- rownames(VIP2D)
    v  <- VIP2D[, comp]
    tim <- sub(paste0("^(.*)", sep, "(.*)$"), "\\2", rn)
    idx_by_time <- split(seq_along(rn), tim)
    top_names <- unlist(lapply(idx_by_time, function(idx) {
      ord <- idx[order(-v[idx])]
      rn[head(ord, min(top_n, length(ord)))]
    }), use.names = FALSE)
    intersect(top_names, rn_target)  # return only hits among our observed Top-N
  }
  
  .permute_top_pvals_parallel <- function(
    df_top, X, Y, comp, top_n, R, ncomp, sep, fit_fun,
    seed = NULL, workers = NULL, cap_blas_threads = FALSE,
    parallel = c("auto","off"), future_packages = NULL, ...
  ) {
    parallel <- match.arg(parallel)
    rn_target <- paste0(df_top$feature, sep, df_top$time)
    counts <- setNames(integer(length(rn_target)), rn_target)
    n <- dim(X)[1]
    
    # cap temporaneo dei thread BLAS/OpenMP 
    if (isTRUE(cap_blas_threads)) {
      if (requireNamespace("withr", quietly = TRUE)) {
        withr::local_envvar(c(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1"))
      } else {
        old_omp <- Sys.getenv("OMP_NUM_THREADS", NA_character_)
        old_mkl <- Sys.getenv("MKL_NUM_THREADS", NA_character_)
        on.exit({
          if (!is.na(old_omp)) Sys.setenv(OMP_NUM_THREADS = old_omp) else Sys.unsetenv("OMP_NUM_THREADS")
          if (!is.na(old_mkl)) Sys.setenv(MKL_NUM_THREADS = old_mkl) else Sys.unsetenv("MKL_NUM_THREADS")
        }, add = TRUE)
        Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
      }
    }
    
    # funzione per una singola permutazione 
    perm_fun <- function(r){
      if (!is.null(seed)) set.seed(seed + r)
      perm <- sample.int(n)
      Yp <- if (length(dim(Y))==3) Y[perm, , , drop = FALSE] else Y[perm, , drop = FALSE]
      VIP2D <- fit_fun(X, Yp, ncomp = ncomp, ...)
      rn <- rownames(VIP2D)
      v  <- VIP2D[, comp]
      tim <- sub(paste0("^(.*)", sep, "(.*)$"), "\\2", rn)
      idx_by_time <- split(seq_along(rn), tim)
      top_names <- unlist(lapply(idx_by_time, function(idx) {
        ord <- idx[order(-v[idx])]
        rn[head(ord, min(top_n, length(ord)))]
      }), use.names = FALSE)
      intersect(top_names, rn_target)
    }
    
    # future.apply
    use_future <- (parallel == "auto") &&
      requireNamespace("future.apply", quietly = TRUE) &&
      requireNamespace("future", quietly = TRUE) &&
      R > 0 && isTRUE(tryCatch(future::nbrOfWorkers() > 1, error = function(e) FALSE))
    
    if (use_future) {
      hits_list <- future.apply::future_lapply(
        seq_len(R), perm_fun,
        future.seed = TRUE,
        future.packages = unique(na.omit(future_packages))  # es: c("mixOmics","<tuoPkg>")
      )
    } else {
      hits_list <- lapply(seq_len(R), perm_fun)
    }
    
    # aggrega
    if (length(hits_list)) {
      all_hits <- unlist(hits_list, use.names = FALSE)
      if (length(all_hits)) {
        tab <- table(all_hits)
        counts[names(tab)] <- as.integer(tab)
      }
    }
    p <- (counts + 1) / (R + 1)
    data.frame(row = names(p), p_perm = as.numeric(p), stringsAsFactors = FALSE)
  }
  
  # inputs & VIP2D
  VIP2D <- if (is.matrix(obj)) obj else if (!is.null(obj$VIP2D)) obj$VIP2D else
    stop("Pass the nplsda_vips object (with $VIP2D) or the VIP2D matrix.")
  if (!is.matrix(VIP2D)) stop("VIP2D must be a matrix.")
  if (is.null(ncomp)) ncomp <- if (!is.matrix(obj) && !is.null(obj$ncomp_used)) obj$ncomp_used else 1
  if (comp < 1 || comp > ncol(VIP2D))
    stop(sprintf("comp=%d out of range (ncol VIP2D = %d)", comp, ncol(VIP2D)))
  
  if (length(dim(X)) != 3) stop("X must be a 3D array (n × p × k).")
  samp_ids   <- dimnames(X)[[1]]
  feat_names <- dimnames(X)[[2]]
  time_names <- dimnames(X)[[3]]
  if (is.null(samp_ids) || is.null(feat_names) || is.null(time_names))
    stop("Please set dimnames(X): [[1]]=samples, [[2]]=features, [[3]]=time.")
  
  groups <- .coerce_groups(groups, sample_ids = samp_ids, case_level = case_level)
  g_ctrl <- levels(groups)[1]; g_case <- levels(groups)[2]
  if (is.null(fit_fun)) fit_fun <- .fit_fun_default
  
  rn <- rownames(VIP2D); if (is.null(rn)) stop("VIP2D must have rownames 'feature_sep_time'.")
  st <- .split_feature_time(rn, sep = sep)
  df <- tibble::tibble(feature = st$feature, time = st$time, VIP = as.numeric(VIP2D[, comp]))
  
#Top-N per timepoint
  df_top <- df |>
    dplyr::group_by(time) |>
    dplyr::slice_max(order_by = VIP, n = top_n, with_ties = FALSE) |>
    dplyr::ungroup()
  
  #  permutation p-values 
  if (permute_R > 0) {
    if (is.null(Y)) stop("Provide Y when permute_R > 0 (needed to permute labels).")
    ptab <- .permute_top_pvals_parallel(
      df_top, X, Y, comp, top_n, permute_R, ncomp, sep, fit_fun,
      seed = seed,
      cap_blas_threads = cap_blas_threads,
      parallel = parallel,
      future_packages = future_packages,
      ...
    )
    df_top <- dplyr::left_join(
      df_top,
      transform(ptab,
                feature = sub(paste0("^(.*)", sep, "(.*)$"), "\\1", row),
                time    = sub(paste0("^(.*)", sep, "(.*)$"), "\\2", row)
      ),
      by = c("feature","time")
    )
    df_top$signif_perm <- !is.na(df_top$p_perm) & df_top$p_perm < alpha_perm
  } else {
    df_top$p_perm <- NA_real_
    df_top$signif_perm <- FALSE
  }
  #  tiles (group means) for selected rows
  rows_list <- lapply(seq_len(nrow(df_top)), function(i) {
    f  <- df_top$feature[i]; tm <- df_top$time[i]
    j  <- match(f, feat_names); tt <- match(tm, time_names)
    if (is.na(j) || is.na(tt)) stop("Feature/time not found in X dimnames: ", f, " / ", tm)
    xvec <- X[, j, tt]
    m_case <- mean(xvec[groups == g_case], na.rm = TRUE)
    m_ctrl <- mean(xvec[groups == g_ctrl], na.rm = TRUE)
    data.frame(
      time = tm, feature = f,
      group = factor(c(g_case, g_ctrl), levels = c(g_case, g_ctrl)),
      value = c(m_case, m_ctrl),
      winner = c(m_case >= m_ctrl, m_ctrl > m_case),
      stringsAsFactors = FALSE
    )
  })
  tiles <- dplyr::bind_rows(rows_list)
  
  # NO-GAPS: per-time y order shared by both panels
  ord <- df_top |>
    dplyr::group_by(time) |>
    dplyr::arrange(dplyr::desc(VIP), .by_group = TRUE) |>
    dplyr::mutate(y_fac = factor(feature, levels = rev(unique(feature))))
  df_top <- dplyr::left_join(df_top, dplyr::select(ord, feature, time, y_fac), by = c("feature","time"))
  tiles  <- dplyr::left_join(tiles,  dplyr::select(ord, feature, time, y_fac), by = c("feature","time"))
  
  #  Left: lollipop (orange bars; green points if significant)
  orange <- "#F28E2B"
  green  <- "#2EAD3A"  # vivid green
  
  p_left <- ggplot2::ggplot(df_top, ggplot2::aes(x = VIP, y = y_fac, group = feature)) +
    ggplot2::geom_segment(ggplot2::aes(x = threshold, xend = VIP, yend = y_fac),
                          linewidth = 0.7, color = orange) +
    ggplot2::geom_point(size = 2.2, color = orange) +
    ggplot2::geom_vline(xintercept = threshold, linetype = 2) +
    ggplot2::facet_wrap(~ time, scales = "free_y") +
    ggplot2::scale_y_discrete(drop = TRUE) +
    ggplot2::labs(title = sprintf("VIP2D — Top %d per timepoint (Component %d)", top_n, comp),
                  x = "VIP score", y = NULL) +
    ggplot2::theme_minimal(base_size = 11)
  
  # overlay: permutation significance -> green point
  if (any(df_top$signif_perm)) {
    p_left <- p_left +
      ggplot2::geom_point(
        data = subset(df_top, signif_perm),
        ggplot2::aes(x = VIP, y = y_fac),
        size = 2.8, color = green
      )
  }
  
  # Right: strip
  if (mode == "winner") {
    p_right <- ggplot2::ggplot(tiles, ggplot2::aes(x = group, y = y_fac)) +
      ggplot2::geom_tile(ggplot2::aes(fill = group, alpha = winner),
                         width = 0.98, height = 0.9, color = NA) +
      ggplot2::scale_fill_manual(values = stats::setNames(c(case_col, control_col), c(g_case, g_ctrl))) +
      ggplot2::scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = loser_alpha), guide = "none") +
      ggplot2::facet_wrap(~ time, scales = "free_y") +
      ggplot2::scale_y_discrete(drop = TRUE) +
      ggplot2::labs(x = NULL, y = NULL, fill = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank())
  } else {
    p_right <- ggplot2::ggplot(tiles, ggplot2::aes(x = group, y = y_fac)) +
      ggplot2::geom_tile(ggplot2::aes(fill = value),
                         width = 0.98, height = 0.9, color = NA) +
      ggplot2::scale_fill_gradient2(low = control_col, mid = "white", high = case_col, midpoint = 0, name = "Mean") +
      ggplot2::facet_wrap(~ time, scales = "free_y") +
      ggplot2::scale_y_discrete(drop = TRUE) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 11) +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank())
  }
  
  #  Combine
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(patchwork::wrap_plots(p_left, p_right, widths = c(4, 1)))
  } else if (requireNamespace("cowplot", quietly = TRUE)) {
    return(cowplot::plot_grid(p_left, p_right, ncol = 2, rel_widths = c(4, 1)))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    return(gridExtra::grid.arrange(p_left, p_right, ncol = 2, widths = c(4, 1)))
  } else {
    stop("Please install one of: patchwork, cowplot, or gridExtra to combine panels.")
  }
}



