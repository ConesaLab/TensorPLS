#' Plot NPLS-DA Scores
#'
#' @param scores_matrix Matrix of scores (from NPLSDAvariates$X)
#' @param expl_var_matrix Matrix with variance explained (with R2.Y column)
#' @param class_vec Factor vector with class labels
#' @param pc1 First component to plot (default = 1)
#' @param pc2 Second component to plot (default = 2)
#' @param main_title Plot title (default = "NPLS-DA Score Plot")
#'
#' @export
plot_nplsda_scores <- function(scores_matrix,
                               expl_var_matrix,
                               class_vec,
                               pc1 = 1,
                               pc2 = 2,
                               variance = FALSE,
                               main_title = "NPLS-DA Score Plot") {
  
  if (pc1 > ncol(scores_matrix) || pc2 > ncol(scores_matrix)) {
    stop(sprintf("Component indices (%d, %d) exceed available components (%d)", 
                 pc1, pc2, ncol(scores_matrix)))
  }
  
  
  # Extract coordinates
  x <- scores_matrix[, pc1]
  y <- scores_matrix[, pc2]
  
  
  if (length(class_vec) != length(x)) {
    stop(sprintf("Length mismatch: class_vec has %d elements, scores have %d observations",
                 length(class_vec), length(x)))
  }
  # Axis labels with R2.Y percentages
  xlab <- sprintf("Comp %d: %.1f%%", pc1, expl_var_matrix[pc1, 3] * 100)
  ylab <- sprintf("Comp %d: %.1f%%", pc2, expl_var_matrix[pc2, 3] * 100)
  
  # Set colors based on classes
  class_levels <- levels(class_vec)
  n_classes <- length(class_levels)
  color_palette <- c("dodgerblue2", "orange2", "green3", "red3", "purple3", "gold2")
  class_colors <- setNames(color_palette[1:n_classes], class_levels)
  cols <- class_colors[as.character(class_vec)]
  
  # Plot with external legend
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(5, 4, 4, 8))
  if(variance){
  plot(x, y,
       pch = 21, 
       bg = cols,
       col = "black",
       xlab = xlab, 
       ylab = ylab,
       main = main_title)
  
  abline(h = 0, v = 0, lty = 2, col = "black")
  }
  else{
    plot(x,y,
         pch = 21,
         bg = cols,
         col = "black",
         xlab = sprintf("Comp %d ", pc1),
         ylab = sprintf("Comp %d ", pc2),
         main = main_title)
    abline (h = 0, v = 0, lty = 2, col = "black")
  }
  # External legend - set xpd=TRUE only for legend
  par(xpd = TRUE)
  legend(x = par("usr")[2] + diff(par("usr")[1:2]) * 0.02,
         y = mean(par("usr")[3:4]),
         legend = names(class_colors),
         pch = 21,
         pt.bg = class_colors,
         col = "black",
         bty = "n",
         title = "Classes")
  
  invisible(NULL)
}
