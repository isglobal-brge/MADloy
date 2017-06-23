#' @method plot LOY
#' @export
plot.LOY <- function(x, labels, colors = c("red", "blue", "darkgreen"), pos.leg = "bottomleft", 
    cex.label = 0.8, print.labels = FALSE, ...) {
    data <- x$data
    if (missing(labels)) {
        if (attr(x, "type") == "Coverage") 
            labels <- rownames(data)
        if (attr(x, "type") == "LRR") 
            labels <- names(data)
    }
    
    
    nclass <- length(unique(x$class))
    
    if (nclass == 1) {
        mycol <- colors[2]
        leg.lab <- c("normal")
        col.lab <- colors[2]
        # alt <- rep(FALSE, length(x$class))
    } else if (nclass == 2) {
        loy <- x$class == "LOY"
        gain <- x$class == "XYY"
        if (sum(loy) > sum(gain)) {
            mycol <- ifelse(loy, colors[1], colors[2])
            leg.lab <- c("LOY", "normal")
            col.lab <- colors[1:2]
            alt <- loy
        } else {
            mycol <- ifelse(gain, colors[3], colors[2])
            leg.lab <- c("normal", "XYY")
            col.lab <- colors[2:3]
            # alt <- gain
        }
    }
    if (nclass == 3) {
        loy <- x$class == "LOY"
        gain <- x$class == "XYY"
        mycol <- ifelse(loy, colors[1], ifelse(gain, colors[3], colors[2]))
        leg.lab <- c("LOY", "normal", "XYY")
        col.lab <- colors
        # alt <- x$class%in%c('LOY', 'XYY')
    }
    
    
    if (attr(x, "type") == "Coverage") {
        ss <- 1:nrow(data)
        xx <- data[, 1]
        yy <- data[, 2]
        d <- xx - yy
        graphics::plot.default(ss, d, type = "n", xlab = "Individuals", ylab = "Mean coverage difference (mY region - Reference)")
        graphics::points(ss, d, col = mycol, pch = 16, ...)
    }
    if (attr(x, "type") == "LRR") {
        ss <- 1:length(data)
        d <- x$data
        graphics::plot.default(ss, d, type = "n", xlab = "Individuals", ylab = "Offset-adjusted trimmed mean normalized mLRR-Y")
        graphics::points(ss, d, col = mycol, pch = 16, ...)
    }
    graphics::legend(pos.leg, leg.lab, pch = 16, col = col.lab, horiz = TRUE, cex = 0.8)
    alt <- x$class %in% c("LOY", "XYY")
    if (any(alt) & print.labels) {
        if (requireNamespace("wordcloud", quietly = TRUE)) {
            wordcloud::textplot(x = ss[alt], y = d[alt], words = tools::file_path_sans_ext(labels[alt]), 
                cex = cex.label, new = FALSE, xlim = c(min(ss), max(ss)), ylim = c(min(d[alt]), 
                  max(d[alt])))
        } else {
            graphics::text(ss[alt], jitter(d[alt]), tools::file_path_sans_ext(labels[alt]), 
                cex = cex.label, adj = 0)
        }
    }
    graphics::abline(h = 0, lty = 2, col = "red")
}

