#' @method plot MADloy
#' @export
plot.MADloy <- function(x, ...) {
  data <- MADloy:::getSummary(x)
  xx <- data[, 1]
  yy <- data[, 2]
  d <- xx - yy
  ss <- 1:nrow(data)
  ref <- GenomeInfoDb::seqnames(x$par$ref.region)
  ref <- ifelse( length(GenomeInfoDb::seqnames(x$par$ref.region)) == 22 , "Autosomes", paste(GenomeInfoDb::seqnames(x$par$ref.region), collapse="_"))
  plot.default(ss, d, type = "n", xlab = "Individuals", 
               ylab = paste0("Mean LRR difference (",GenomeInfoDb::seqnames(x$par$target.region), " - ", ref, ")"), ...)
  points(ss, d, pch = 16)
  abline(h=0, lty=2, col="red")
} 
