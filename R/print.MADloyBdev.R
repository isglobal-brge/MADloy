#' @export
print.MADloyBdev <- function(x, ...) {
  cat("Object of class MADloyBdev \n")
  cat("---------------------------- \n")
  cat("Number of processed samples:", length(x$par$files), "\n")
  cat("Human Genome build:", resBdev$par$hg, "\n")
  cat("Top BAF Threshold:", resBdev$par$top, "\n")
  cat("Bot BAF Threshold:", resBdev$par$bot, "\n")
  cat("\n")
  # TODO adapt to MADloy processed objects
} 