#' Get median values from MADseqLOY or MADloy objects
#' 
#' This function retrieves the median values from the MADloy or MADseqLOY
#' objects and prepares the data to be analyzed for LOY events.
#' 
#' @param object A MADloy or MADseqLOY object from the \code{MADloy} or 
#'   \code{MADseqLOY} functions.
#'   
#' @return A data table with the results for all samples in columns
#' @examples
#' \dontrun{
#' getMedian(resMADloy)
#' getMedian(resMADseqLOY)}
getMedian <- function(object) {
  
  if (inherits(object, "MADseqLOY")) {
    targetAvg <- sapply(object$target, "[[", "medianTargetCoverage")
    refAvg.1 <- sapply(object$reference.1, "[[", "medianTargetCoverage")
    refAvg.2 <- sapply(object$reference.2, "[[", "medianTargetCoverage")
    avg <- cbind(targetAvg, refAvg.1, refAvg.2)
    targetChr <- as.character(GenomeInfoDb::seqnames(object$par$target.region))
    refChr.1 <- as.character(GenomeInfoDb::seqnames(object$par$ref.region.1))
    refChr.2 <- as.character(GenomeInfoDb::seqnames(object$par$ref.region.2))
    colnames(avg) <- c(paste0("medianCov_", targetChr), paste0("medianCov_", refChr.1), paste0("medianCov_", refChr.2))
  }
  
  if (inherits(object, "MADloy")) {
    targetAvg <- sapply(object$target, "[[", "median")
    refAvg <- sapply(object$reference, "[[", "median")
    avg <- cbind(targetAvg, refAvg)
    targetChr <- as.character(GenomeInfoDb::seqnames(object$par$target.region))
    refChr <- as.character(GenomeInfoDb::seqnames(object$par$ref.region))
    if ( length(refChr) > 1 ) refChr <- paste(refChr, collapse = "_")
    colnames(avg) <- c(paste0("medianLRR_", targetChr), paste0("medianLRR_", refChr))
  }
  return(avg)
} 
