#' @export
print.MADloy <- function(x, ...) {
  cat("Object of class MADloy \n")
  cat("---------------------------- \n")
  cat("Number of processed samples:", length(x$par$files), "\n")
  cat("Target region:", paste0("chr", as.character(GenomeInfoDb::seqnames(x$par$target.region)), ":", BiocGenerics::start(x$par$target.region), "-", BiocGenerics::end(x$par$target.region)), "\n")
  if (identical(as.character(seqnames(x$par$ref.region)), as.character(1:22))){
    cat("Reference region(s): Autosomal chromosomes\n")
  } else {
    cat("Reference region(s):", paste0("chr", as.character(GenomeInfoDb::seqnames(x$par$ref.region)), ":", BiocGenerics::start(x$par$ref.region), "-", BiocGenerics::end(x$par$ref.region)), "\n")
  }
  cat("Offset (median LRR value in msY):", round(x$par$offset,2), "\n")
  cat("\n")
  
} 
