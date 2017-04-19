#' Process MAD file to obtain Bdev
#' 
#' Process a MAD file to retrieve the B deviation in a region of
#' interest.
#' 
#' @param file File path of the MAD file to be processed.
#' @param query Region of interest to check the LRR in rangedData format.
#' @param rsCol Column of the MAD file with the name of the SNP.
#' @param ChrCol Column of the MAD file with the chromosome information.
#' @param PosCol Column of the MAD file with the position information.
#' @param BAFCol Column of the MAD file with the B Allele Frequency information.
#' @param top Superior treshold to consider an heterozygous allele.
#' @param bot Inferior treshold to consider an heterozygous allele.
#' @return A list with the Bdev for the query value.
#' @import data.table
#' @examples
#' \dontrun{
#' processBdevMAD(file=file.path, query=rangedDataROI, rsCol=1, chrCol=2, PosCol=3, LRRCol=4)}
processBdevMAD <- function(file, query, rsCol, ChrCol, PosCol, BAFCol, top, bot) {
  
  dat <- data.table::fread(file, showProgress = FALSE, sep = "\t")
  if (!all(c("Name", "Chr", "Position", "B.Allele.Freq") %in% colnames(dat))) {
    stop("Wrong MAD file header in ", file)
  }
  data.table::setnames(dat, colnames(dat[, c(rsCol, ChrCol, PosCol, BAFCol), with = F]), 
    c("Name", "Chr", "Position", "B.Allele.Freq"))
  Bdevsummary <- list()
  sel <- dat[ which(dat$Chr == as.character(GenomeInfoDb::seqnames(query[1])) & dat$Position > BiocGenerics::start(query[1]) & dat$Position < BiocGenerics::end(query[1]) & dat$B.Allele.Freq <= top & dat$B.Allele.Freq >= bot)]
  Bdevsummary$PAR1$Bdev <- mean(abs(0.5 - sel$B.Allele.Freq), na.rm = T)
  Bdevsummary$PAR1$Bdevsd <- sd(abs(0.5 - sel$B.Allele.Freq), na.rm = T)
  Bdevsummary$PAR1$LRR <- mean(abs(0.5 - sel$Log.R.Ratio), na.rm = T)
  Bdevsummary$PAR1$LRRsd <- sd(abs(0.5 - sel$Log.R.Ratio), na.rm = T)
  Bdevsummary$PAR1$n <- nrow(sel)
  sel <- dat[ which(dat$Chr == as.character(GenomeInfoDb::seqnames(query[2])) & dat$Position > BiocGenerics::start(query[2]) & dat$Position < BiocGenerics::end(query[2]) & dat$B.Allele.Freq <= top & dat$B.Allele.Freq >= bot)]
  Bdevsummary$PAR2$Bdev <- mean(abs(0.5 - sel$B.Allele.Freq), na.rm = T)
  Bdevsummary$PAR2$Bdevsd <- sd(abs(0.5 - sel$B.Allele.Freq), na.rm = T)
  Bdevsummary$PAR2$LRR <- mean(abs(0.5 - sel$Log.R.Ratio), na.rm = T)
  Bdevsummary$PAR2$LRRsd <- sd(abs(0.5 - sel$Log.R.Ratio), na.rm = T)
  Bdevsummary$PAR2$n <- nrow(sel)
  sel <- dat[ which(dat$Chr == as.character(GenomeInfoDb::seqnames(query[3])) & as.numeric(dat$Position) > BiocGenerics::start(query[3]) & as.numeric(dat$Position) < BiocGenerics::end(query[3]))]
  Bdevsummary$p$Bdev <- mean(abs(0.5 - sel$B.Allele.Freq), na.rm = T)
  Bdevsummary$p$Bdevsd <- sd(abs(0.5 - sel$B.Allele.Freq), na.rm = T)
  Bdevsummary$p$LRR <- mean(abs(0.5 - sel$Log.R.Ratio), na.rm = T)
  Bdevsummary$p$LRRsd <- sd(abs(0.5 - sel$Log.R.Ratio), na.rm = T)
  Bdevsummary$p$n <- nrow(sel)
  sel <- dat[ which(dat$Chr == as.character(GenomeInfoDb::seqnames(query[4])) & as.numeric(dat$Position) > BiocGenerics::start(query[4]) & as.numeric(dat$Position) < BiocGenerics::end(query[4]))]
  Bdevsummary$q$Bdev <- mean(abs(0.5 - sel$B.Allele.Freq), na.rm = T)
  Bdevsummary$q$Bdevsd <- sd(abs(0.5 - sel$B.Allele.Freq), na.rm = T)
  Bdevsummary$q$LRR <- mean(abs(0.5 - sel$Log.R.Ratio), na.rm = T)
  Bdevsummary$q$LRRsd <- sd(abs(0.5 - sel$Log.R.Ratio), na.rm = T)
  Bdevsummary$q$n <- nrow(sel)
  
  return(Bdevsummary)
} 
