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
#' @return A list with the mean of the LRR for the query value.
#' @import data.table
#' @examples
#' \dontrun{
#' processMAD(file=file.path, query=rangedDataROI, rsCol=2, chrCol=3, PosCol=4, LRRCol=5)}
processMAD <- function(file, query, rsCol, ChrCol, PosCol, LRRCol) {
  
  dat <- data.table::fread(file, showProgress = FALSE, sep = "\t")
  if (!all(c("Name", "Chr", "Position", "Log.R.Ratio") %in% colnames(dat))) {
    stop("Wrong MAD file header in ", file)
  }
  data.table::setnames(dat, colnames(dat[, c(rsCol, ChrCol, PosCol, LRRCol), with = F]), 
    c("Name", "Chr", "Position", "Log.R.Ratio"))
  LRRmedian <- list()
  if( length(query) >1 ) {
    LRRmedian$median <- median(dat$Log.R.Ratio[which(dat$Chr == as.character(GenomeInfoDb::seqnames(query)) & dat$Position > BiocGenerics::start(query) & dat$Position < BiocGenerics::end(query))], na.rm = T)
  } else {
    LRRmedian$median <- median(dat$Log.R.Ratio[which(apply(sapply(query, function(x) dat$Chr == as.character(GenomeInfoDb::seqnames(x)) & dat$Position > BiocGenerics::start(x) & dat$Position < BiocGenerics::end(x)), 1, any))], na.rm = T)
  }
  
  return(LRRmedian)
} 
