#' @method plot MADseqLOY
#' @export
plot.MADseqLOY <- function(x, ...) {
  par(mfrow=c(1,2))
  data <- MADloy:::getSummary(x)
  xx <- data[, 1]
  yy <- data[, 2]
  d <- xx - yy
  ss <- 1:nrow(data)
  plot.default(ss, d, type = "n", xlab = "Individuals", ylab = paste0("Mean coverage difference (",GenomeInfoDb::seqnames(x$par$target.region), " - ", GenomeInfoDb::seqnames(x$par$ref.region), ")"), ...)
  points(ss, d, pch = 16)
} 
