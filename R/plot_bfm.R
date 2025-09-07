#' Pareto frontier (min complexity, max Fit%)
#'
#' Plot all models (points) and the Pareto frontier (step line) where
#' no model can improve Fit% without increasing total complexity.
#'
#' @param tab A data frame with columns \code{SumComp} (P+Q+R) and \code{Fitper}
#'   plus the usual \code{P,Q,R}. Typically \code{out$summary} from
#'   \code{bestfittedmodel()}.
#'
#' @return A \code{ggplot} object (invisible).
#' @examples
#' \donttest{
#' # out <- bestfittedmodel(X, 2, 2, 4)
#' # plot_pareto_front(out$summary)
#' }
#' @export
plot_pareto_front <- function(tab,plotName = NULL) {
  stopifnot(all(c("SumComp", "Fitper") %in% names(tab)))
  df <- tab[order(tab$SumComp, -tab$Fitper), , drop = FALSE]
  
  # keep = running "upper envelope" across increasing SumComp
  keep <- logical(nrow(df))
  best <- -Inf
  for (i in seq_len(nrow(df))) {
    fi <- df$Fitper[i]
    if (is.na(fi)) next
    if (fi >= best) {
      keep[i] <- TRUE
      best <- fi
    }
  }
  
  g <- ggplot2::ggplot(df, ggplot2::aes(SumComp, Fitper)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_step(
      data = df[keep, , drop = FALSE],
      ggplot2::aes(y = Fitper),
      direction = "vh"
    ) +
    ggplot2::geom_point(
      data = df[keep, , drop = FALSE],
      size = 2
    ) +
    ggplot2::labs(
      title = plotName,
      x = "P + Q + R",
      y = "Fit (%)"
    ) +
    ggplot2::theme_minimal()
  
  print(g)
  invisible(g)
}

#' Fit% heatmap (facet by one mode)
#'
#' Heatmap of explained variance (Fit%) over a 2D grid of components while
#' faceting by the third mode. Draws a bold outline around the best cell
#' in each facet (ties broken by lower \code{SumComp}).
#'
#' @param tab A data frame with columns \code{P,Q,R,Fitper,SumComp}.
#' @param facet Which mode to facet by: \code{"R"}, \code{"P"}, or \code{"Q"}.
#' @param show_values Logical; print Fit% inside tiles (1 decimal).
#' @param outline_best Logical; outline the best tile in each facet.
#'
#' @return A \code{ggplot} object (invisible).
#' @examples
#' \donttest{
#' # out <- bestfittedmodel(X, 2, 2, 4)
#' # plot_heatmap_components(out$summary, facet = "R")
#' }
#' @importFrom rlang .data
#' @export

plot_heatmap_components <- function(tab, facet = "R",
                                    show_values = TRUE,
                                    outline_best = TRUE
                                    ) {
  stopifnot(all(c("P","Q","R","Fitper","SumComp") %in% names(tab)))
  df  <- subset(tab, !is.na(Fitper))
  rng <- range(df$Fitper, na.rm = TRUE)
  
  # choose which mode to facet, and which go on x/y
  if (facet == "R") {
    x <- "P"; y <- "Q"; f <- "R"
    xlab <- "P (mode-1 components)"; ylab <- "Q (mode-2 components)"
  } else if (facet == "P") {
    x <- "Q"; y <- "R"; f <- "P"
    xlab <- "Q (mode-2 components)"; ylab <- "R (mode-3 components)"
  } else if (facet == "Q") {
    x <- "P"; y <- "R"; f <- "Q"
    xlab <- "P (mode-1 components)"; ylab <- "R (mode-3 components)"
  } else stop("facet must be 'R', 'P' or 'Q'")
  
  # best row per facet (tie-break: higher Fit%, then lower SumComp)
  if (outline_best) {
    ord  <- df[order(df[[f]], -df$Fitper, df$SumComp), ]
    best <- ord[!duplicated(ord[[f]]), ]
  }
  
  library(ggplot2)
  
  g <- ggplot(df, aes(x = factor(.data[[x]]),
                      y = factor(.data[[y]]),
                      fill = Fitper)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_viridis_c(name = "Fit%", limits = rng,
                         breaks = scales::pretty_breaks(5),
                         na.value = "grey90") +
    facet_wrap(as.formula(paste("~", f)), nrow = 1, labeller = label_both) +
    labs(title = "Fit% over P×Q grid with fixed facet",
         x = xlab, y = ylab) +
    coord_fixed() +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          strip.text = element_text(face = "bold"))
  
  if (show_values) {
    g <- g + geom_text(aes(label = sprintf("%.1f", Fitper)), size = 3)
  }
  
  if (outline_best) {
    # draw a bold outline around the best tile
    g <- g + geom_tile(data = best,
                       aes(x = factor(.data[[x]]), y = factor(.data[[y]])),
                       inherit.aes = FALSE, fill = NA, color = "black",
                       linewidth = 1.2)
  }
  
  g
}


