#' Plots the SNP data of chromosome X for an individual
#' 
#' This function draws the LRR and BAF data of chromosome X for an individual analyzed with madloy or checkBdev
#' 
#' @param x A MADloy object from the \code{MADloy} functions.
#' @param sample The identifier of a sample in the MADloy object.
#' @param rsCol Column of the MAD file with the name of the SNP.
#' @param ChrCol Column of the MAD file with the chromosome information.
#' @param PosCol Column of the MAD file with the position information.
#' @param LRRCol Column of the MAD file with the LRR information.
#' @param BAFCol Column of the MAD file with the BAF information.
#' @param offset offset value to adjust msY region LRR. By default is set to 0.
#' @param ... Any other graphical parameter.
#' 
#' @return A data table with the results for all samples in columns
#' @export
#' @examples
#' \dontrun{
#' plotindSNP(resMADloy, "SAMPLE")}
plotIndSNPX <- function(x, sample, rsCol=1, ChrCol=2, PosCol=3, LRRCol=4, BAFCol=5, offset=0, ...) {
    lrr2ploidy <- function(x) 2 * exp(3 * x/2)
    ploidy2lrr <- function(x) 2 * log(x/2)/3
    if (inherits(x, "MADloy") | inherits(x, "MADloyBdev")) {
    samples <- x$Bdev$par$files
    paths <- x$Bdev$par$path
    regions <- x$Bdev$par$regions
    if(is.numeric(sample)){
      ss <- samples[sample]
      pp <- paths[sample]
    } else {
      if (sample%in%tools::file_path_sans_ext(samples)) {
        ii <- which(sample==tools::file_path_sans_ext(samples))
        ss <- samples[ii]
        pp <- paths[ii]
      } else stop("The selected sample has not been processed")
    }
    dat <- data.frame(fread(file.path(pp, ss), header = TRUE))
    o <- grep("^X", (dat[, ChrCol]))
    xy <- any(c("XY", "25") %in% unique(dat$Chr))
    if (xy){
      xysel <- c(grep("^XY", (dat[, ChrCol])), grep("^23", (dat[, ChrCol])))
      datxy <- dat[xysel, ]
      selPAR1 <- datxy[, PosCol] <= as.numeric(regions[regions$chromosome == "X" & regions$type == "PAR1"]$end)
      selPAR2 <- datxy[, PosCol] >= as.numeric(regions[regions$chromosome == "X" & regions$type == "PAR2"]$start) & datxy[, PosCol] <= as.numeric(regions[regions$chromosome == "X" & regions$type == "PAR2"]$end)
    }
    dat <- dat[o, ]
    ## Scale LRR of the msY region to have ploidy = 2 = 1(Y) + 1(constant)
    graphics::par(mar = c(5, 5, 4, 5) + 0.1)
    graphics::plot(dat[, PosCol], dat[, LRRCol]+offset, ylim = c(-2, 2), las = 1, pch =".", cex = 2, col = 1, ylab = "",xlab = "", main = "", xaxt="n", yaxt="n", xlim=c(1, as.numeric(regions[regions$chromosome == "X" & regions$type == "PAR2"]$end)), ...)
    if (xy){
      graphics::points(datxy[selPAR1, PosCol], datxy[selPAR1, LRRCol], pch =".", cex = 2, col=1)
      graphics::points(datxy[selPAR2, PosCol], datxy[selPAR2, LRRCol], pch =".", cex = 2, col=1)
    }
    graphics::axis( 2, at = c(-2:2), labels = c("-2.0", "-1.0", "0.0", "1.0", "2.0"), las = 1, col = "black", col.axis = "black")
    #abline(h=x$par$offset, col="orange", lty=2, lwd=2)
    graphics::par(new = TRUE)
    graphics::plot(dat[, PosCol], dat[, BAFCol], col = 2, pch = ".", cex = 2, ylab = "", xlab = "", main = "", axes = F, xlim=c(1, as.numeric(regions[regions$chromosome == "X" & regions$type == "PAR2"]$end)), ...)
    if (xy){
      graphics::points(datxy[selPAR1, PosCol], datxy[selPAR1, BAFCol], pch =".", cex = 2, col=2)
      graphics::points(datxy[selPAR2, PosCol], datxy[selPAR2, BAFCol], pch =".", cex = 2, col=2)
    }
    graphics::abline(h=0.5, col=8)
    graphics::abline(h=c(0.33, 0.66), col=8)
    graphics::mtext("LRR", side = 2, col = "black", line = 2.5, adj = 0.5)
    graphics::mtext("BAF", side = 4, col = "red", line = 2.5, adj = 0.5)
    graphics::axis( 4, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.0", 0.25, 0.5, 0.75, "1.0"), las = 1, col = "black", col.axis = "red")
    xaxis <- seq(0, max(dat[, PosCol], na.rm=T), by=25000000)
    graphics::axis( 1, at = xaxis, labels = as.character(round(xaxis/1000000, 1)), las = 1, col = "black", col.axis = "black")
    graphics::mtext("position (Mb)", side = 1, col = "black", line = 2.5)
    graphics::rect(regions[regions$chromosome == "X" & regions$type == "PAR1"]$start, -0.01, regions[regions$chromosome == "X" & regions$type == "PAR1"]$end, -0.03, col="Blue", border=NA)
    graphics::rect(regions[regions$chromosome == "X" & regions$type == "PAR1"]$start, 1.01, regions[regions$chromosome == "X" & regions$type == "PAR1"]$end, 1.03, col="Blue", border=NA)
    graphics::abline(v=c(regions[regions$chromosome == "X" & regions$type == "PAR1"]$start, regions[regions$chromosome == "X" & regions$type == "PAR1"]$end), lty=2, col="blue")
    graphics::rect(regions[regions$chromosome == "X" & regions$type == "PAR2"]$start, -0.01, regions[regions$chromosome == "X" & regions$type == "PAR2"]$end, -0.03, col="Blue", border=NA)
    graphics::rect(regions[regions$chromosome == "X" & regions$type == "PAR2"]$start, 1.01, regions[regions$chromosome == "X" & regions$type == "PAR2"]$end, 1.03, col="Blue", border=NA)
    graphics::abline(v=c(regions[regions$chromosome == "X" & regions$type == "PAR2"]$start, regions[regions$chromosome == "X" & regions$type == "PAR2"]$end), lty=2, col="blue")
    graphics::rect(regions[regions$chromosome == "X" & regions$type == "XTR"]$start, -0.01, regions[regions$chromosome == "X" & regions$type == "XTR"]$end, -0.03, col="Orange", border=NA)
    graphics::rect(regions[regions$chromosome == "X" & regions$type == "XTR"]$start, 1.01, regions[regions$chromosome == "X" & regions$type == "XTR"]$end, 1.03, col="Orange", border=NA)
    graphics::abline(v=c(regions[regions$chromosome == "X" & regions$type == "XTR"]$start, regions[regions$chromosome == "X" & regions$type == "XTR"]$end), lty=2, col="Orange")
    graphics::title(sample)
    graphics::title("Chromosome X", line = 0.3)
  } else {
    stop("The object supplied is not a MADloy or MADloyBdev object.")
  }
}
