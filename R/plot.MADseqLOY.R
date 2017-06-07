#' @method plot MADseqLOY
#' @export
plot.MADseqLOY <- function(x, ...) {
  graphics::par(mfrow=c(1,2))
  data <- getSummary(x)
  xx <- data[, 1]
  yy <- data[, 2]
  d <- xx - yy
  ss <- 1:nrow(data)
  graphics::plot.default(ss, d, type = "n", xlab = "Individuals", ylab = paste0("Mean coverage difference (",GenomeInfoDb::seqnames(x$par$target.region), " - ", GenomeInfoDb::seqnames(x$par$ref.region), ")"), ...)
  graphics::points(ss, d, pch = 16)
} 
