#' Get summary values from MADseqLOY or MADloy objects
#' 
#' This function retrieves the summary values from the MADloy or MADseqLOY
#' objects and prepares the data to be analyzed for LOY events.
#' 
#' @param object A MADloy or MADseqLOY object from the \code{MADloy} or 
#'   \code{MADseqLOY} functions.
#'   
#' @return A data table with the results for all samples in columns
#' @examples
#' \dontrun{
#' getSummary(resMADloy)
#' getSummary(resMADseqLOY)}
getSummary <- function(object) {
    
    if (inherits(object, "MADseqLOY")) {
        targetAvg <- sapply(object$target, "[[", "medianTargetCoverage")
        refAvg <- sapply(object$reference, "[[", "medianTargetCoverage")
        avg <- cbind(targetAvg, refAvg)
        targetChr <- as.character(GenomeInfoDb::seqnames(object$par$target.region))
        refChr <- as.character(GenomeInfoDb::seqnames(object$par$ref.region))
        colnames(avg) <- c(paste0("summaryCov_", targetChr), paste0("summaryCov_", 
            refChr))
    }
    
    if (inherits(object, "MADloy")) {
        targetAvg <- sapply(object$target, "[[", "summary")
        refAvg <- sapply(object$reference, "[[", "summary")
        avg <- cbind(targetAvg, refAvg)
        targetChr <- as.character(GenomeInfoDb::seqnames(object$par$target.region))
        refChr <- as.character(GenomeInfoDb::seqnames(object$par$ref.region))
        if (identical(refChr, as.character(1:22))) 
            refChr <- "Autosomes"
        colnames(avg) <- c(paste0("summaryLRR_", targetChr), paste0("summaryLRR_", 
            refChr))
    }
    return(avg)
}
