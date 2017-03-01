#' @export
print.MADseqLOY <- function(x, ...) {
  cat("\n")
  cat("Object of class MADseqLOY \n")
  cat("---------------------------- \n")
  cat("Number of processed samples:", length(x$par$files), "\n")
  cat("Target region:", paste0("chr", as.character(GenomeInfoDb::seqnames(x$par$target.region)), ":", BiocGenerics::start(x$par$target.region), "-", BiocGenerics::end(x$par$target.region)), "\n")
  cat("Reference region 1:", paste0("chr", as.character(GenomeInfoDb::seqnames(x$par$ref.region)), ":", BiocGenerics::start(x$par$ref.region), "-", BiocGenerics::end(x$par$ref.region)), "\n")
  cat("\n")
  # TODO adapt to MADseqloy processed objects
} 