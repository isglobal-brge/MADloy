#' @method plot MADseqLOY
#' @export
plot.MADseqLOY <- function(x, ...) {
  par(mfrow=c(1,2))
  data <- MADloy:::getMedian(x)
  xx <- data[, 1]
  yy.1 <- data[, 2]
  d.1 <- xx - yy.1
  ss <- 1:nrow(data)
  plot.default(ss, d.1, type = "n", xlab = "Individuals", ylab = paste0("Median coverage difference (",GenomeInfoDb::seqnames(x$par$target.region), " - ", GenomeInfoDb::seqnames(x$par$ref.region), ")"), ...)
  points(ss, d.1, pch = 16)
} 
