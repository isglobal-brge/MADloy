#' @method plot MADloy
#' @export
plot.MADloy <- function(x, ...) {
par(mfrow=c(1,2))
  data <- MADloy:::getMedian(x)
  xx <- data[, 1]
  yy.1 <- data[, 2]
  yy.2 <- data[, 3]
  d.1 <- xx - yy.1
  d.2 <- xx - yy.2
  ss <- 1:nrow(data)
  plot.default(ss, d.1, type = "n", xlab = "Individuals", ylab = paste0("Median LRR difference (",GenomeInfoDb::seqnames(x$par$target.region), " - ", GenomeInfoDb::seqnames(x$par$ref.region.1), ")"), ...)
  points(ss, d.1, pch = 16)
  plot.default(ss, d.2, type = "n", xlab = "Individuals", ylab = paste0("Median LRR difference (",GenomeInfoDb::seqnames(x$par$target.region), " - ", GenomeInfoDb::seqnames(x$par$ref.region.2), ")"), ...)
  points(ss, d.2, pch = 16)
} 
