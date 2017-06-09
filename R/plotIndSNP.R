#' Plots the SNP data of an individual
#' 
#' This function draws the LRR and BAF data for an individual
#' 
#' @param x A MADloy object from the \code{MADloy} functions or A single file path (APT platform and MAD platform).
#' @param sample The identifier of a sample in the MADloy object when x is a MADloy object.
#' @param rsCol Column of the MAD file with the name of the SNP.
#' @param ChrCol Column of the MAD file with the chromosome information.
#' @param PosCol Column of the MAD file with the position information.
#' @param LRRCol Column of the MAD file with the LRR information.
#' @param BAFCol Column of the MAD file with the BAF information.
#' @param offset offset applied to the msY region.
#' @param hg Human genome build version. It can be "hg18", "hg19" or "GRCh38". Set by default to "hg18". Only used when a file path is given.
#' @param ... Any other graphical parameter.
#' 
#' @return A data table with the results for all samples in columns
#' @export
#' @examples
#' \dontrun{
#' plotindSNP(resMADloy, "SAMPLE")}
plotIndSNP <- function(x, sample, rsCol=1, ChrCol=2, PosCol=3, LRRCol=4, BAFCol=5, offset=0, hg="hg18", ...) {
  if (inherits(x, "MADloy") | inherits(x, "MADloyBdev")) {
    
    samples <- x$par$files
    paths <- x$par$path
    regions <- x$par$regions
    if(is.numeric(sample)){
      ss <- samples[sample]
      pp <- paths[sample]
    } else {
      if (sample%in%tools::file_path_sans_ext(samples)) {
        ii <- grep(sample, samples)
        ss <- samples[ii]
        pp <- paths[ii]
      } else stop("The selected sample has not been processed")
    }
    dat <- data.frame(fread(file.path(pp, ss), header = TRUE))
    o <- grep("^Y", (dat[, ChrCol]))
    xy <- any(c("XY", "25") %in% unique(dat$Chr))
    if (xy){
      xysel <- c(grep("^XY", (dat[, ChrCol])), grep("^23", (dat[, ChrCol])))
      datxy <- dat[xysel, ]
    }
    dat <- dat[o, ]
    graphics::par(mar = c(5, 5, 4, 5) + 0.1)
    graphics::plot.default(dat[, PosCol], dat[, LRRCol]+offset, ylim = c(-5, 5), las = 1, pch =".", cex = 2, col = 1, ylab = "",xlab = "", main = "", xaxt="n", yaxt="n", xlim=c(1, as.numeric(regions[3, 4])), ...)
    if (xy){
      graphics::points(datxy[datxy[, PosCol] <= as.numeric(regions[2, 4]), PosCol], datxy[datxy[, PosCol] <= as.numeric(regions[2, 4]), LRRCol], pch =".", cex = 2, col=1)
      graphics::points(datxy[datxy[, PosCol] >= as.numeric(regions[4, 3]) & datxy[, PosCol] <= as.numeric(regions[4, 4]), PosCol]-as.numeric(regions[3,4])+as.numeric(regions[3,3]), datxy[datxy[, PosCol] >= as.numeric(regions[4, 3]) & datxy[, PosCol] <= as.numeric(regions[4, 4]), LRRCol], pch =".", cex = 2, col=1)
    }
    graphics::axis( 2, at = c(-5:5), labels = c("-5.0", "-4.0", "-3.0", "-2.0", "-1.0", "0.0", "1.0", "2.0", "3.0", "4.0", "5.0"), las = 1, col = "black", col.axis = "black")
    graphics::par(new = TRUE)
    graphics::plot.default(dat[, PosCol], dat[, BAFCol], col = 2, pch = ".", cex = 2, ylab = "", xlab = "", main = "", axes = F, xlim=c(1, as.numeric(regions[3, 4])), ...)
    if (xy){
      graphics::points(datxy[datxy[, PosCol] <= as.numeric(regions[2, 4]), PosCol], datxy[datxy[, PosCol] <= as.numeric(regions[2, 4]), BAFCol], pch =".", cex = 2, col=2)
      graphics::points(datxy[datxy[, PosCol] >= as.numeric(regions[4, 3]) & datxy[, PosCol] <= as.numeric(regions[4, 4]), PosCol]-as.numeric(regions[3,4])+as.numeric(regions[3,3]), datxy[datxy[, PosCol] >= as.numeric(regions[4, 3]) & datxy[, PosCol] <= as.numeric(regions[4, 4]), BAFCol], pch =".", cex = 2, col=1)
    }
    graphics::abline(h=0.5, col=8)
    graphics::abline(h=c(0.33, 0.66), col=8)
    graphics::mtext("LRR", side = 2, col = "black", line = 2.5, adj = 0.5)
    graphics::mtext("BAF", side = 4, col = "red", line = 2.5, adj = 0.5)
    graphics::axis( 4, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.0", 0.25, 0.5, 0.75, "1.0"), las = 1, col = "black", col.axis = "red")
    xaxis <- seq(0, max(dat[, PosCol], na.rm=T), by=25000000)
    graphics::axis( 1, at = xaxis, labels = as.character(round(xaxis/1000000, 1)), las = 1, col = "black", col.axis = "black")
    graphics::mtext("position (Mb)", side = 1, col = "black", line = 2.5)
    graphics::rect(6671498, -0.01, 22919969, -0.03, col="green", border=NA)
    graphics::rect(6671498, 1.01, 22919969, 1.03, col="green", border=NA)
    graphics::abline(v=c(6671498, 22919969), lty=2, col="green")
    graphics::rect(1, -0.01, 2709520, -0.03, col="Blue", border=NA)
    graphics::rect(1, 1.01, 2709520, 1.03, col="Blue", border=NA)
    graphics::abline(v=c(1, 2709520), lty=2, col="blue")
    graphics::rect(57443438, -0.01, 57772954, -0.03, col="Blue", border=NA)
    graphics::rect(57443438, 1.01, 57772954, 1.03, col="Blue", border=NA)
    graphics::abline(v=c(57443438, 57772954), lty=2, col="blue")
    graphics::title(sample)
    graphics::title("Chromosome Y", line = 0.3)
  } else {
    regions <- fread(system.file("extdata", "references", paste0(hg, ".par.regions"), package = "MADloy"), header=T, skip=1, colClasses = c("character", "character", "numeric", "numeric"), showProgress = FALSE)
    if(missing(sample)) sample <- tools::file_path_sans_ext(basename(x))
    dat <- data.frame(fread(x, header = TRUE))
    o <- grep("^Y", (dat[, ChrCol]))
    xy <- any(c("XY", "25") %in% unique(dat$Chr))
    if (xy){
      xysel <- c(grep("^XY", (dat[, ChrCol])), grep("^23", (dat[, ChrCol])))
      datxy <- dat[xysel, ]
    }
    dat <- dat[o, ]
    graphics::par(mar = c(5, 5, 4, 5) + 0.1)
    graphics::plot.default(dat[, PosCol], dat[, LRRCol]+offset, ylim = c(-5, 5), las = 1, pch =".", cex = 2, col = 1, ylab = "",xlab = "", main = "", xaxt="n", yaxt="n", xlim=c(1, as.numeric(regions[3, 4])), ...)
    if (xy){
      graphics::points(datxy[datxy[, PosCol] <= as.numeric(regions[2, 4]), PosCol], datxy[datxy[, PosCol] <= as.numeric(regions[2, 4]), LRRCol], pch =".", cex = 2, col=1)
      graphics::points(datxy[datxy[, PosCol] >= as.numeric(regions[4, 3]) & datxy[, PosCol] <= as.numeric(regions[4, 4]), PosCol]-as.numeric(regions[3,4])+as.numeric(regions[3,3]), datxy[datxy[, PosCol] >= as.numeric(regions[4, 3]) & datxy[, PosCol] <= as.numeric(regions[4, 4]), LRRCol], pch =".", cex = 2, col=1)
    }
    graphics::axis( 2, at = c(-5:5), labels = c("-5.0", "-4.0", "-3.0", "-2.0", "-1.0", "0.0", "1.0", "2.0", "3.0", "4.0", "5.0"), las = 1, col = "black", col.axis = "black")
    graphics::par(new = TRUE)
    graphics::plot.default(dat[, PosCol], dat[, BAFCol], col = 2, pch = ".", cex = 2, ylab = "", xlab = "", main = "", axes = F, xlim=c(1, as.numeric(regions[3, 4])), ...)
    if (xy){
      graphics::points(datxy[datxy[, PosCol] <= as.numeric(regions[2, 4]), PosCol], datxy[datxy[, PosCol] <= as.numeric(regions[2, 4]), BAFCol], pch =".", cex = 2, col=2)
      graphics::points(datxy[datxy[, PosCol] >= as.numeric(regions[4, 3]) & datxy[, PosCol] <= as.numeric(regions[4, 4]), PosCol]-as.numeric(regions[3,4])+as.numeric(regions[3,3]), datxy[datxy[, PosCol] >= as.numeric(regions[4, 3]) & datxy[, PosCol] <= as.numeric(regions[4, 4]), BAFCol], pch =".", cex = 2, col=1)
    }
    graphics::abline(h=0.5, col=8)
    graphics::abline(h=c(0.33, 0.66), col=8)
    graphics::mtext("LRR", side = 2, col = "black", line = 2.5, adj = 0.5)
    graphics::mtext("BAF", side = 4, col = "red", line = 2.5, adj = 0.5)
    graphics::axis( 4, at = c(0, 0.25, 0.5, 0.75, 1), labels = c("0.0", 0.25, 0.5, 0.75, "1.0"), las = 1, col = "black", col.axis = "red")
    xaxis <- seq(0, max(dat[, PosCol], na.rm=T), by=25000000)
    graphics::axis( 1, at = xaxis, labels = as.character(round(xaxis/1000000, 1)), las = 1, col = "black", col.axis = "black")
    graphics::mtext("position (Mb)", side = 1, col = "black", line = 2.5)
    graphics::rect(6671498, -0.01, 22919969, -0.03, col="green", border=NA)
    graphics::rect(6671498, 1.01, 22919969, 1.03, col="green", border=NA)
    graphics::abline(v=c(6671498, 22919969), lty=2, col="green")
    graphics::rect(1, -0.01, 2709520, -0.03, col="Blue", border=NA)
    graphics::rect(1, 1.01, 2709520, 1.03, col="Blue", border=NA)
    graphics::abline(v=c(1, 2709520), lty=2, col="blue")
    graphics::rect(57443438, -0.01, 57772954, -0.03, col="Blue", border=NA)
    graphics::rect(57443438, 1.01, 57772954, 1.03, col="Blue", border=NA)
    graphics::abline(v=c(57443438, 57772954), lty=2, col="blue")
    graphics::title(sample)
    graphics::title("Chromosome Y", line = 0.3)
  }
}
