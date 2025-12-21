#' Plot Block X Mode 3
#'
#' @param factors  Object returned from compute_npls_factors
#'                 ($FactorsX$Mode3)
#' @param pcs      vector of integer containing the components to plot (default c(1,2))
#' @param edge     vecctors containing size of border (default(c(xlimSize,ylimSize)))
#' @return         Invisibile, dataframe of points plotted 
#' @export
plot_nplsda_blockX_mode3 <- function(factors, pcs = c(1, 2),edge = c(0.2,0.3)) {
  stopifnot(is.list(factors), !is.null(factors$FactorsX$Mode3))
  M3 <- factors$FactorsX$Mode3
  stopifnot(ncol(M3) >= max(pcs))
  
  lab <- rownames(M3)
  if (is.null(lab)) lab <- paste0("t", seq_len(nrow(M3)))
  
  x <- M3[, pcs[1]]
  y <- M3[, pcs[2]]
  
  xlim <- c(min(x), max(x) + abs(edge[1] * max(x)))
  ylim <- c(min(y), max(y) + abs(edge[2] * max(y)))
  
  df <- data.frame(time = lab, x = x, y = y, stringsAsFactors = FALSE)
  
  if (requireNamespace("ggplot2", quietly = TRUE) &&
      requireNamespace("ggrepel", quietly = TRUE)) {
    
    gg <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point(shape = 21, size = 3,
                          colour = "black", fill = "dodgerblue2") +
      ggrepel::geom_text_repel(ggplot2::aes(label = time),
                               colour = "dodgerblue3", size = 4,
                               min.segment.length = 0, seed = 1) +
      ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
      ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
      ggplot2::scale_x_continuous(limits = xlim,
                                  name = paste("Component", pcs[1])) +
      ggplot2::scale_y_continuous(limits = ylim,
                                  name = paste("Component", pcs[2])) +
      ggplot2::ggtitle("Block X Mode 3: NPLS-DA") +
      ggplot2::theme_minimal(base_size = 12)
    print(gg)
    
  } else {
    # Fallback base R con jitter per ridurre le sovrapposizioni
    oldpar <- par(no.readonly = TRUE); on.exit(par(oldpar))
    par(mar = c(5, 5, 3, 2))
    plot(x, y, type = "p", col = "black", pch = 21, bg = "dodgerblue2",
         xlab = paste("Component", pcs[1]),
         ylab = paste("Component", pcs[2]),
         main = "Block X Mode 3: NPLS-DA",
         xlim = xlim, ylim = ylim)
    abline(v = 0, h = 0, lty = 2, col = "grey40")
    # piccola spinta radiale per evitare collisioni
    ang  <- atan2(y - mean(y), x - mean(x))
    text(x + 0.02 * diff(range(x)) * cos(ang),
         y + 0.02 * diff(range(y)) * sin(ang),
         labels = lab, col = "dodgerblue3", cex = 0.95, pos = 3)
  }
  
  invisible(df)
}
#fullarrayunionPos
#plot_nplsda_blockX_mode3(PosFactors,edge = c(0.2,0.3))
#plot_nplsda_blockX_mode3(NegFactors,edge = c(0.2,0.3))
#plot_nplsda_blockX_mode3(GEFactors,edge = c(0.2,0.3))
#plot_nplsda_blockX_mode3(GCTOFFactors,edge = c(0.2,0.3))
#PosFactors = compute_npls_factors(X = fullarrayunionPos,Y = outcomedummyarray136,ncomp = 2 )
#NegFactors = compute_npls_factors(X = Fullarraynplsda_vipsIntersectionNeg2DModel2, Y = outcomedummyarray136,ncomp = 2) 
#GEFactors = compute_npls_factors(X = fullarrayintersectionGElist, Y = outcomedummyarray136,ncomp = 2 )
#GCTOFFactors = compute_npls_factors(X = fullarrayintersectionGCTOFX2DMode2l, Y = outcomedummyarray136, ncomp = 2)
#MetabolomicsFactors = compute_npls_factors(X = Allmetabolomics136,outcomedummyarray136,ncomp = 2)
#plot_nplsda_blockX_mode3(MetabolomicsFactors,edge = c(0.2,0.3))
#MetabolomicsFactors
#PosFactors
#GEFactors = compute_npls_factors(X = fullarrayintersectionGElist, Y = outcomedummyarray136,ncomp = 2 )
#plot_nplsda_blockX_mode3(GEFactors,edge = c(0.2,0.3))
