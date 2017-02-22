#' @method plot MADloy
#' @export
plot.MADloy <- function(x, ...) {
  data <- MADloy:::getMedian(x)
  xx <- data[, 1]
  yy.1 <- data[, 2]
  d.1 <- xx - yy.1
  ss <- 1:nrow(data)
  ref <- GenomeInfoDb::seqnames(x$par$ref.region)
  ref <- ifelse( length(GenomeInfoDb::seqnames(x$par$ref.region)) == 22 , "Autosomes", paste(GenomeInfoDb::seqnames(x$par$ref.region), collapse="_"))
  plot.default(ss, d.1, type = "n", xlab = "Individuals", ylab = paste0("Median LRR difference (",GenomeInfoDb::seqnames(x$par$target.region), " - ", ref, ")"), ...)
  points(ss, d.1, pch = 16)
} 
