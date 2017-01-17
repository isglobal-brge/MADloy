
plotIndSNP <- function(x, sample, ...) {

  samples <- x$par$files
  paths <- x$par$path
  if(is.numeric(sample)){
    ss <- samples[sample]
    pp <- paths[sample]
  }  
  else {
   if (sample%in%tools::file_path_sans_ext(samples)) {
    ii <- grep(sample, samples)
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
  
  data.table::setnames(dat, colnames(dat[, c(rsCol, ChrCol, PosCol, LRRCol), with = F]), 
                       c("Name", "Chr", "Position", "Log.R.Ratio"))
  queryA <- unlist(strsplit(x = "chrY:7000000-30000000", split = "[:, -]", perl = T))
  subsetA <- GenomicRanges::GRanges(seqnames = gsub("chr", "", queryA[1]), 
                                    ranges = IRanges::IRanges(start = as.numeric(queryA[2]), end = as.numeric(queryA[3])))
  
  lrr.target<- dat$Log.R.Ratio[which(dat$Chr == as.character(GenomeInfoDb::seqnames(subsetA)) & dat$Position > 
                          BiocGenerics::start(subsetA) & dat$Position < BiocGenerics::end(subsetA))]
  pos.target <- dat$Position[which(dat$Chr == as.character(GenomeInfoDb::seqnames(subsetA)) & dat$Position > 
                                        BiocGenerics::start(subsetA) & dat$Position < BiocGenerics::end(subsetA))]
  
  
  plot(pos.target, lrr.target, ylab="LRR", xlab="Position (Mb) - Chr Y", type="n", ...)
  rect(7e6, max(lrr.target, na.rm=TRUE), 25e6, min(lrr.target, na.rm=TRUE), col="MistyRose", border="white")
  points(pos.target, lrr.target, pch=16, cex=0.7, col="blue")
  title(tools::file_path_sans_ext(ss))
  abline(h=0, col="red", lty=2, lwd=3)
}