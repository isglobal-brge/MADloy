#' Plots the SNP LRR data of an individual
#' 
#' This function draws the LRR data for an individual
#' 
#' @param x A MADloy object from the \code{MADloy} functions.
#' @param sample The identifier of a sample in the MADloy object.
#' @param rsCol Column of the MAD file with the name of the SNP.
#' @param ChrCol Column of the MAD file with the chromosome information.
#' @param PosCol Column of the MAD file with the position information.
#' @param LRRCol Column of the MAD file with the LRR information.
#' @param ... Other parameters to be passed to the graphics::plot function.   
#' @return A data table with the results for all samples in columns
#' @export
#' @examples
#' \dontrun{
#' plotindLRR(resMADloy, "SAMPLE")}
plotIndLRR <- function(x, sample, rsCol=1, ChrCol=2, PosCol=3, LRRCol=4, ...) {
  if (inherits(x, "MADloy")) {
   samples <- x$par$files
   paths <- x$par$path
   
   if(is.numeric(sample)){
    ss <- samples[sample]
    pp <- paths[sample]
   }  
   else {
    if (sample%in%tools::file_path_sans_ext(samples)) {
     ii <- which(sample==tools::file_path_sans_ext(samples))
     ss <- samples[ii]
     pp <- paths[ii]
    }  
    else
     stop("The selected sample has not been processed")
    }
  
  dat <- fread(file.path(pp, ss), header = TRUE)
  rsCol <- x$par$cols[1]
  ChrCol <- x$par$cols[2]
  PosCol <- x$par$cols[3]
  LRRCol <- x$par$cols[4]
  tt <- tools::file_path_sans_ext(ss)
  }
  
  else if (is.data.frame(x)) {
    x.i <- x[, c(rsCol, ChrCol, PosCol, LRRCol)]
    if (!is.numeric(x.i[,3])) {
      x.i[,3] <- as.numeric(as.character(x.i[,3]))
      warning("\n Genomic position information has changed to numeric values. Please, check \n")
    }  
    dat <- as.data.table(x.i)
    tt <- names(x.i)[4]
    rsCol <- 1
    ChrCol <- 2
    PosCol <- 3
    LRRCol <- 4
  }
  
    
  else {
    dat <- fread(x, header=TRUE)
    tt <- tools::file_path_sans_ext(basename(x))
  }  
    
  data.table::setnames(dat, colnames(dat[, c(rsCol, ChrCol, PosCol, LRRCol), with = F]), 
                       c("Name", "Chr", "Position", "Log.R.Ratio"))
  queryA <- unlist(strsplit(x = "chrY:2694521-59034049", split = "[:, -]", perl = T))
  subsetA <- GenomicRanges::GRanges(seqnames = gsub("chr", "", queryA[1]), 
                                    ranges = IRanges::IRanges(start = as.numeric(queryA[2]), end = as.numeric(queryA[3])))
  
  lrr.target<- dat$Log.R.Ratio[which(dat$Chr == as.character(GenomeInfoDb::seqnames(subsetA)) & dat$Position > 
                          BiocGenerics::start(subsetA) & dat$Position < BiocGenerics::end(subsetA))]
  pos.target <- dat$Position[which(dat$Chr == as.character(GenomeInfoDb::seqnames(subsetA)) & dat$Position > 
                                        BiocGenerics::start(subsetA) & dat$Position < BiocGenerics::end(subsetA))]
  
  if (max(lrr.target, na.rm=TRUE)<1){
    graphics::plot(pos.target, lrr.target, ylab="LRR", xlab="Position (Mb) - Chr Y", type="n", 
                   ylim=c(min(lrr.target, na.rm=TRUE), 1), ...)
    uu <- graphics::par("usr")
  }  
  else {
    graphics::plot(pos.target, lrr.target, ylab="LRR", xlab="Position (Mb) - Chr Y", type="n", ...)
    uu <- graphics::par("usr")
  }
  graphics::rect(7e6, uu[4], 25e6, uu[3], col="MistyRose")
  graphics::points(pos.target, lrr.target, pch=16, cex=0.7, col="brown")
  graphics::title(tt)
  graphics::abline(h=0, col="red", lty=2, lwd=2)
  graphics::abline(h=x$par$offset, col="orange", lty=2, lwd=2)
  hh <- x$target[[ss]]$summary
  graphics::segments(7e6, hh, 25e6, col="blue", lwd=3)
}
