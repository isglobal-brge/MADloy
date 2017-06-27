#' Process MAD file
#' 
#' Process a MAD file to retrieve the median log R ratio of a region of
#' interest.
#' 
#' @param file File path of the MAD file to be processed.
#' @param query Region of interest to check the LRR in rangedData format.
#' @param rsCol Column of the MAD file with the name of the SNP.
#' @param ChrCol Column of the MAD file with the chromosome information.
#' @param PosCol Column of the MAD file with the position information.
#' @param LRRCol Column of the MAD file with the LRR information.
#' @param trim the fraction (0 to 0.5) of probes to be trimmed when summaryzing LRR. 
#' This argument tries to control the effect of having CNVs across genome. Default is 0.1.
#' @return A list with the mean of the LRR for the query value.
#' @import data.table
#' @examples
#' \dontrun{
#' processMAD(file=file.path, query=rangedDataROI, rsCol=2, chrCol=3, PosCol=4, LRRCol=5)}
processMAD <- function(file, query, rsCol, ChrCol, PosCol, LRRCol, trim=0.1) {
  
  dat <- data.table::fread(file, showProgress = FALSE, sep = "\t")
  
  data.table::setnames(dat, colnames(dat[, c(rsCol, ChrCol, PosCol, LRRCol), with = F]), 
    c("Name", "Chr", "Position", "Log.R.Ratio"))
  LRRsummary <- list()
  if( length(query) == 1 ) {
    LRRsummary$summary <- mean(dat$Log.R.Ratio[which(dat$Chr == as.character(GenomeInfoDb::seqnames(query)) & dat$Position > BiocGenerics::start(query) & dat$Position < BiocGenerics::end(query))], na.rm = T, trim=trim)
    LRRsummary$sd <- sd(dat$Log.R.Ratio[which(dat$Chr == as.character(GenomeInfoDb::seqnames(query)) & dat$Position > BiocGenerics::start(query) & dat$Position < BiocGenerics::end(query))], na.rm = T)
  } else if( identical(as.character(GenomeInfoDb::seqnames(query)), as.character(1:22))){
    LRRsummary$summary <- mean(dat$Log.R.Ratio[which(dat$Chr %in% as.character(1:22))], na.rm = T, trim=trim)
    LRRsummary$sd <- sd(dat$Log.R.Ratio[which(dat$Chr %in% as.character(1:22))], na.rm = T)
  } else {
    LRRsummary$summary <- mean(dat$Log.R.Ratio[which(apply(sapply(query, function(x) dat$Chr == as.character(GenomeInfoDb::seqnames(x)) & dat$Position > BiocGenerics::start(x) & dat$Position < BiocGenerics::end(x)), 1, any))], na.rm = T, trim=trim)
    LRRsummary$sd <- sd(dat$Log.R.Ratio[which(apply(sapply(query, function(x) dat$Chr == as.character(GenomeInfoDb::seqnames(x)) & dat$Position > BiocGenerics::start(x) & dat$Position < BiocGenerics::end(x)), 1, any))], na.rm = T)
  }
  
  return(LRRsummary)
} 
