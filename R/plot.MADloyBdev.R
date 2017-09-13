#' @method plot MADloyBdev
#' @export
plot.MADloyBdev <- function(x, ...) {
    sel <- x$class$balanced != "balancedPAR" & !is.na(x$class$balanced)
    filled <- rep(1, sum(sel))
    filled[(x$class$balanced != "balancedPAR" & !is.na(x$class$balanced))[sel]] <- 20
    
    gplots::plotCI(x = unlist(data.frame(t(sapply(x$Bdev, "[[", "PAR1")))$Pl), 
        uiw = unlist(data.frame(t(sapply(x$Bdev, "[[", "PAR1")))$Plsd), lty = 2, 
        xaxt = "n", pch = filled, gap = 0, ylim = c(0, 3), minbar = 0, maxbar = 3, 
        main = "Estimated ploidy trimmed mean Bdeviation in PAR1 and PAR2 regions", ylab = "Estimated ploidy", 
        xlab = "getLOY classification")
    
    graphics::par(new = TRUE)
    
    gplots::plotCI(x = unlist(data.frame(t(sapply(x$Bdev, "[[", "PAR2")))$Pl), 
        uiw = unlist(data.frame(t(sapply(x$Bdev, "[[", "PAR2")))$Plsd), lty = 2, 
        xaxt = "n", pch = filled, gap = 0, ylim = c(0, 3), minbar = 0, maxbar = 3, 
        col = "red", ylab = "", xlab = "")
    
    abline(h=2, lty=2, col="orange")
    
    graphics::axis(1, at = c(1:length(x$Bdev)), labels = FALSE)
    
    graphics::text(x = c(1:length(x$Bdev)), y = graphics::par()$usr[3] - 0.125, labels = x$class$orig, 
                   srt = 45, adj = 1, xpd = TRUE, cex = 0.75)
    
    graphics::par(new = TRUE)
    
    gplots::plotCI(x = unlist(data.frame(t(sapply(x$Bdev, "[[", "PAR1")))$Bdev),
        xaxt="n", pch = 20, gap = 0, ylim=c(0, 0.5), minbar=0, maxbar=3, ylab ="",
        xlab ="", yaxt="n")

    graphics::par(new = TRUE)
    
    gplots::plotCI(x = unlist(data.frame(t(sapply(x$Bdev, "[[", "PAR2")))$Bdev),
                   xaxt = "n", pch = 20, gap = 0, ylim = c(0, 0.5), minbar = 0,
                   maxbar = 3, ylab = "", xlab = "", yaxt = "n", col = "red")
    
    graphics::axis(4, seq(0, 0.5, 0.1))
    graphics::mtext("Bdeviation", side = 4, line = 2.5, adj = 0.5)
    
    abline(h=0.06, lty=2, col="gray")
    
    graphics::legend("topright", legend = c("Ploidy PAR1", "Ploidy PAR2", "Bdev PAR1", "Bdev PAR2"), pch = c(21, 21, 20, 20), col = c("black", 
        "red", "black", "red"))
    
    
}
