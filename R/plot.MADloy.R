#' Plot of MADloy object
#' 
#' This function plots the summarized LRR of a set of samples
#' 
#' @param x a MADloy object.
#' @param labels a character vector with the samples labels.
#' @param print.labels should sample labels be showed.
#' @param threshold threshold to draw sample labels.
#' @param cex.label size of the labels.
#' @param ... Other graphical parameters.
#' @method plot MADloy
#' @export
plot.MADloy <- function(x, labels, print.labels=FALSE, 
                        threshold = -0.6, cex.label=0.8, ...) {
  
  mLRRY <- x$mLRRY$mLRRY
  ss <- 1:length(mLRRY)
  d <- mLRRY
  if (missing(labels))
    labels <- names(d)
  ref <- GenomeInfoDb::seqnames(x$par$ref.region)
  ref <- ifelse( length(GenomeInfoDb::seqnames(x$par$ref.region)) == 22 , "Autosomes", paste(GenomeInfoDb::seqnames(x$par$ref.region), collapse="_"))
  graphics::plot.default(ss, d, type = "n", xlab = "Individuals", 
               ylab = "Trimmed mean normalized mLRR-Y", ...)
  graphics::points(ss, d, pch = 16)
  graphics::abline(h=0, lty=2, col="red")
 
  alt <- d <= threshold
  alt[is.na(alt)] <- FALSE
  
  if (any(alt) & print.labels) {
    if (requireNamespace("wordcloud", quietly = TRUE) & sum(alt)>1) {
     wordcloud::textplot(x = ss[alt], y = d[alt], 
                         words = tools::file_path_sans_ext(labels[alt]), 
                         cex = cex.label, 
                         new = FALSE, xlim=c(min(ss), max(ss)), 
                         ylim=c(min(d[alt]), max(d[alt])))
   } else {
    graphics::text(ss[alt], jitter(d[alt]), tools::file_path_sans_ext(labels[alt]),
         cex = cex.label, adj = 0)
   }
  }
} 
