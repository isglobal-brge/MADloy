#' Get summary values from a MADloy object
#' 
#' This function retrieves the summary values from a MADloy 
#' object and prepares the data to be analyzed for LOY events.
#' 
#' @param object A MADloy object from the \code{MADloy} function.
#'   
#' @return A data table with the results for all samples in columns.
#' @examples
#' \dontrun{
#' getSummary(resMADloy)}
getSummary <- function(object) {

  if (inherits(object, "MADloy")) {
    targetAvg <- sapply(object$target, "[[", "summary")
    refAvg <- sapply(object$reference, "[[", "summary")
    avg <- cbind(targetAvg, refAvg)
    targetChr <- as.character(GenomeInfoDb::seqnames(object$par$target.region))
    refChr <- as.character(GenomeInfoDb::seqnames(object$par$ref.region))
    if ( identical(sort(refChr), sort(as.character(1:22))) ) refChr <- "Autosomes"
    else if ( length(refChr) > 1 )refChr <- "Ref_Region"
    colnames(avg) <- c(paste0("summaryLRR_", targetChr), paste0("summaryLRR_", refChr))
  }
  return(avg)
} 
