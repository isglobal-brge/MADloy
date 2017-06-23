#' @method plot checkSex
#' @export
plot.checkSex <- function(x, ...) {
    
    lrr2ploidy <- function(x) 2 * exp(3 * x/2)
    ploidy2lrr <- function(x) 2 * log(x/2)/3
    
    tmp <- x$data
    tmp$X <- ploidy2lrr(lrr2ploidy(tmp$X) + (2 - lrr2ploidy(x$par$offsetX)))
    tmp$Y <- tmp$Y - x$par$offsetY
    
    ## Calculate LRR - Ploidy lines
    ploidy <- c(0.001, 1, 2, 3, 4, 5)
    lrrploidy <- ploidy2lrr(ploidy)
    posMeans <- apply(cbind(lrrploidy[1:5], lrrploidy[2:6]), 1, mean)
    
    ## Draw the plot
    graphics::plot(tmp, col = x$class, xlim = c(lrrploidy[1], lrrploidy[5]), ylim = c(lrrploidy[1], 
        lrrploidy[5]), ylab = "trimmed mean chrY LRR", xlab = "trimmed mean chrX LRR", 
        pch = 20)
    
    graphics::abline(v = lrrploidy[2:5], lty = 2)
    graphics::abline(h = lrrploidy[2:5], lty = 2)
    graphics::text(x = lrrploidy[2:5], y = lrrploidy[1] + 0.1, labels = c("X", "XX", "XXX", 
        "XXXX"), srt = 90, pos = 4, offset = 0.7)
    graphics::text(x = lrrploidy[1] + 0.1, y = lrrploidy[2:5], labels = c("Y", "YY", "YYY", 
        "YYYY"), pos = 3, offset = 0.1)
    graphics::legend("bottomleft", legend = c("Females", "Males"), pch = 20, col = c(1, 2))
}
