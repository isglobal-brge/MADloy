#' @export
print.MADloyBdev <- function(x, ...) {
  cat("Object of class MADloyBdev \n")
  cat("---------------------------- \n")
  cat("Number of processed samples:", length(x$par$files), "\n")
  cat("Human Genome build:", x$par$hg, "\n")
  cat("Top BAF Threshold:", x$par$top, "\n")
  cat("Bot BAF Threshold:", x$par$bot, "\n")
  cat("\n")
  # TODO adapt to MADloy processed objects
} 