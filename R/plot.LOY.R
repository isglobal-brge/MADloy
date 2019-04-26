#' @method plot LOY
#' @export
plot.LOY <- function(x, labels, colors = c("red", "blue", "darkgreen"), pos.leg = "bottomleft", 
    cex.label = 0.8, print.labels = FALSE, method="MADloy", ...) {
  
  if (method == "MADloy") {
    class <- x$res$MADloy
  } else if (method == "Fosberg") {
    class <- x$res$Fosberg
  } else stop("method should be either 'MADloy' or 'Fosberg'")
  
  data <- x$data
  
  if (missing(labels)) {
    labels <- tools::file_path_sans_ext(rownames(x$res))
  }
    
  nclass <- length(levels(class))
    
    if (nclass == 1) {
        mycol <- colors[2]
        leg.lab <- c("normal")
        col.lab <- colors[2]
        # alt <- rep(FALSE, length(x$class))
    } else if (nclass == 2) {
        loy <- class == "LOY"
        gain <- class == "XYY"
        if (sum(loy, na.rm=T) > sum(gain, na.rm=T)) {
            mycol <- ifelse(loy, colors[1], colors[2])
            leg.lab <- c("LOY", "normal")
            col.lab <- colors[1:2]
            alt <- loy
        } else {
            mycol <- ifelse(gain, colors[3], colors[2])
            leg.lab <- c("normal", "XYY")
            col.lab <- colors[2:3]
        }
    } else if (nclass == 3) {
        loy <- class == "LOY"
        gain <- class == "XYY"
        mycol <- ifelse(loy, colors[1], ifelse(gain, colors[3], colors[2]))
        leg.lab <- c("LOY", "normal", "XYY")
        col.lab <- colors
    }
    
    if (attr(data, "type") == "LRR") {
        ss <- 1:length(data)
        d <- x$data
        graphics::plot.default(ss, d, type = "n", xlab = "Individuals", 
                               ylab = "Offset-adjusted trimmed mean normalized mLRR-Y",
                               ...)
        graphics::points(ss, d, col = mycol, pch = 16, ...)
    }
  
    graphics::legend(pos.leg, leg.lab, pch = 16, col = col.lab, horiz = TRUE, cex = 0.8)
    alt <- class %in% c("LOY", "XYY")
    if (any(alt) & print.labels) {
        if (requireNamespace("wordcloud", quietly = TRUE)) {
            wordcloud::textplot(x = ss[alt], y = d[alt], 
                                words = tools::file_path_sans_ext(labels[alt]), 
                                cex = cex.label, new = FALSE, 
                                xlim = c(min(ss), max(ss)), ylim = c(min(d[alt]),
                                                                     max(d[alt])))
        } else {
            graphics::text(ss[alt], jitter(d[alt]), 
                           tools::file_path_sans_ext(labels[alt]),
                           cex = cex.label, adj = 0)
        }
    }
    graphics::abline(h = 0, lty = 2, col = "red")
}

