#' Process MAD file to obtain Bdev
#' 
#' Process a MAD file to retrieve the B deviation in a region of
#' interest.
#' 
#' @param file File path of the MAD file to be processed.
#' @param regions Regions of chromosomes X and Y in data.frame. See system.file("extdata", "references") for more info on these regions.
#' @param rsCol Column of the MAD file with the name of the SNP.
#' @param ChrCol Column of the MAD file with the chromosome information.
#' @param PosCol Column of the MAD file with the position information.
#' @param LRRCol Column of the MAD file with the Log R Ratio information.
#' @param BAFCol Column of the MAD file with the B Allele Frequency information.
#' @param top Superior treshold to consider an heterozygous allele.
#' @param bot Inferior treshold to consider an heterozygous allele.
#' @param trim the fraction (0 to 0.5) of probes to be trimmed when summaryzing LRR. 
#' This argument tries to control the effect of having CNVs in the analyzed region. Default is 0.1.
#' @return A list with the Bdev for the query value.
#' @import data.table
#' @examples
#' \dontrun{
#' processBdevMAD(file=file.path, query=rangedDataROI, rsCol=1, chrCol=2, PosCol=3, LRRCol=4)}
processBdevMAD <- function(file, regions, rsCol, ChrCol, PosCol, LRRCol, BAFCol, top, 
    bot, trim) {
    lrr2ploidy <- function(x) 2 * exp(3 * x/2)
    ploidy2lrr <- function(x) 2 * log(x/2)/3
    dat <- data.table::fread(file, showProgress = FALSE, sep = "\t")
    xy <- any(c("XY", "25") %in% unique(dat$Chr))
    data.table::setnames(dat, colnames(dat[, c(rsCol, ChrCol, PosCol, LRRCol, BAFCol), 
        with = F]), c("Name", "Chr", "Position", "Log.R.Ratio", "B.Allele.Freq"))
    Bdevsummary <- list()
    if (xy){
      Bdevsummary$check <- "PAR"
      # PAR1
      sel <- dat[which(dat$Chr %in% c("XY", "25") & dat$Position > regions[regions$chromosome == "X" & 
              regions$type == "PAR1"]$start & dat$Position < regions[regions$chromosome == "X" & 
              regions$type == "PAR1"]$end)]
      selhet <- sel$B.Allele.Freq >= bot & sel$B.Allele.Freq <= top & !is.na(sel$B.Allele.Freq)
      Bdevsummary$PAR1$Bdev <- mean(abs(0.5 - sel$B.Allele.Freq[selhet]), na.rm = T)
      Bdevsummary$PAR1$Bdevsd <- stats::sd(abs(0.5 - sel$B.Allele.Freq[selhet]), na.rm = T)
      Bdevsummary$PAR1$Pl <- mean(lrr2ploidy(sel$Log.R.Ratio), na.rm = T, trim = trim)
      Bdevsummary$PAR1$Plsd <- stats::sd(lrr2ploidy(sel$Log.R.Ratio), na.rm = T)
      Bdevsummary$PAR1$n <- sum(selhet)
      Bdevsummary$PAR1$N <- nrow(sel)
      # PAR2
      sel <- dat[which(dat$Chr %in% c("XY", "25") & 
                         dat$Position > regions[regions$chromosome == "X" & regions$type == "PAR2"]$start & 
                         dat$Position < regions[regions$chromosome == "X" & regions$type == "PAR2"]$end)]
      if(nrow(sel) == 0) { 
          stop("No data found in PAR2 region") }
      selhet <- sel$B.Allele.Freq >= bot & sel$B.Allele.Freq <= top & !is.na(sel$B.Allele.Freq)
      Bdevsummary$PAR2$Bdev <- mean(abs(0.5 - sel$B.Allele.Freq[selhet]), na.rm = T)
      Bdevsummary$PAR2$Bdevsd <- stats::sd(abs(0.5 - sel$B.Allele.Freq[selhet]), na.rm = T)
      Bdevsummary$PAR2$Pl <- mean(lrr2ploidy(sel$Log.R.Ratio), na.rm = T, trim = trim)
      Bdevsummary$PAR2$Plsd <- stats::sd(lrr2ploidy(sel$Log.R.Ratio), na.rm = T)
      Bdevsummary$PAR2$n <- sum(selhet)
      Bdevsummary$PAR2$N <- nrow(sel)
      # X
      selPAR1 <- dat$Chr == "X" &
                       dat$Position > regions[regions$chromosome == "X" & regions$type == "PAR1"]$start & 
                       dat$Position < regions[regions$chromosome == "X" & regions$type == "PAR1"]$end &
                       !is.na(dat$Chr) & !is.na(dat$Position)
      selPAR2 <- dat$Chr == "X" &
                       dat$Position > regions[regions$chromosome == "X" & regions$type == "PAR2"]$start & 
                       dat$Position < regions[regions$chromosome == "X" & regions$type == "PAR2"]$end &
                       !is.na(dat$Chr) & !is.na(dat$Position)
      selXTR <- dat$Chr == "X" &
                      dat$Position > regions[regions$chromosome == "X" & regions$type == "XTR"]$start & 
                      dat$Position < regions[regions$chromosome == "X" & regions$type == "XTR"]$end &
                      !is.na(dat$Chr) & !is.na(dat$Position)
      sel <- dat[dat$Chr == "X" & !(selPAR1 | selPAR2 | selXTR) & !is.na(dat$Chr) & !is.na(dat$Position)]
      selhet <- sel$B.Allele.Freq >= bot & sel$B.Allele.Freq <= top & !is.na(sel$B.Allele.Freq)
      Bdevsummary$fsX$Bdev <- mean(abs(0.5 - sel$B.Allele.Freq[selhet]), na.rm = T)
      Bdevsummary$fsX$Bdevsd <- stats::sd(abs(0.5 - sel$B.Allele.Freq[selhet]), na.rm = T)
      Bdevsummary$fsX$Pl <- mean(lrr2ploidy(sel$Log.R.Ratio), na.rm = T, trim = trim)
      Bdevsummary$fsX$Plsd <- stats::sd(lrr2ploidy(sel$Log.R.Ratio), na.rm = T)
      Bdevsummary$fsX$n <- sum(selhet)
      Bdevsummary$fsX$N <- nrow(sel)
      # XTR
      sel <- dat[which(dat$Chr %in% c("X", "23") & 
                         dat$Position > regions[regions$chromosome == "X" & regions$type == "XTR"]$start & 
                         dat$Position < regions[regions$chromosome == "X" & regions$type == "XTR"]$end)]
      selhet <- sel$B.Allele.Freq >= bot & sel$B.Allele.Freq <= top & !is.na(sel$B.Allele.Freq)
      Bdevsummary$XTR$Bdev <- mean(abs(0.5 - sel$B.Allele.Freq[selhet]), na.rm = T)
      Bdevsummary$XTR$Bdevsd <- stats::sd(abs(0.5 - sel$B.Allele.Freq[selhet]), na.rm = T)
      Bdevsummary$XTR$Pl <- mean(lrr2ploidy(sel$Log.R.Ratio), na.rm = T, trim = trim)
      Bdevsummary$XTR$Plsd <- stats::sd(lrr2ploidy(sel$Log.R.Ratio), na.rm = T)
      Bdevsummary$XTR$n <- sum(selhet, na.rm=T)
      Bdevsummary$XTR$N <- nrow(sel)  
    } else {
      warning("No PAR1 or PAR2 probes are available in file", file, " to measure Bdeviation. Using only XTR region to measure Bdeviation.")
      Bdevsummary$check <- "XTR"
      # X
      selPAR1 <- dat$Chr == "X" &
        dat$Position > regions[regions$chromosome == "X" & regions$type == "PAR1"]$start & 
        dat$Position < regions[regions$chromosome == "X" & regions$type == "PAR1"]$end &
        !is.na(dat$Chr) & !is.na(dat$Position)
      selPAR2 <- dat$Chr == "X" &
        dat$Position > regions[regions$chromosome == "X" & regions$type == "PAR2"]$start & 
        dat$Position < regions[regions$chromosome == "X" & regions$type == "PAR2"]$end &
        !is.na(dat$Chr) & !is.na(dat$Position)
      selXTR <- dat$Chr == "X" &
        dat$Position > regions[regions$chromosome == "X" & regions$type == "XTR"]$start & 
        dat$Position < regions[regions$chromosome == "X" & regions$type == "XTR"]$end &
        !is.na(dat$Chr) & !is.na(dat$Position)
      sel <- dat[dat$Chr == "X" & !(selPAR1 | selPAR2 | selXTR) & !is.na(dat$Chr) & !is.na(dat$Position)]
      selhet <- sel$B.Allele.Freq >= bot & sel$B.Allele.Freq <= top & !is.na(sel$B.Allele.Freq)
      Bdevsummary$fsX$Bdev <- mean(abs(0.5 - sel$B.Allele.Freq[selhet]), na.rm = T)
      Bdevsummary$fsX$Bdevsd <- stats::sd(abs(0.5 - sel$B.Allele.Freq[selhet]), na.rm = T)
      Bdevsummary$fsX$Pl <- mean(lrr2ploidy(sel$Log.R.Ratio), na.rm = T, trim = trim)
      Bdevsummary$fsX$Plsd <- stats::sd(lrr2ploidy(sel$Log.R.Ratio), na.rm = T)
      Bdevsummary$fsX$n <- sum(selhet)
      Bdevsummary$fsX$N <- nrow(sel)
      # XTR
      sel <- dat[which(dat$Chr %in% c("X", "23") & 
                         dat$Position > regions[regions$chromosome == "X" & regions$type == "XTR"]$start & 
                         dat$Position < regions[regions$chromosome == "X" & regions$type == "XTR"]$end)]
      selhet <- sel$B.Allele.Freq >= bot & sel$B.Allele.Freq <= top
      Bdevsummary$XTR$Bdev <- mean(abs(0.5 - sel$B.Allele.Freq[selhet]), na.rm = T)
      Bdevsummary$XTR$Bdevsd <- stats::sd(abs(0.5 - sel$B.Allele.Freq[selhet]), na.rm = T)
      Bdevsummary$XTR$Pl <- mean(lrr2ploidy(sel$Log.R.Ratio), na.rm = T, trim = trim)
      Bdevsummary$XTR$Plsd <- stats::sd(lrr2ploidy(sel$Log.R.Ratio), na.rm = T)
      Bdevsummary$XTR$n <- sum(selhet, na.rm=T)
      Bdevsummary$XTR$N <- nrow(sel)  
    }
    return(Bdevsummary)
}
