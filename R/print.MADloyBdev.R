#' @export
print.MADloyBdev <- function(x, ...) {
  cat("Object of class MADloyBdev \n")
  cat("---------------------------- \n")
  cat("Number of balanced pq samples:", sum(x$class$adjusted_p == "balancedpq"), "\n")
  cat("Number of unbalanced pq samples:", sum(x$class$adjusted_p != "balancedpq"), "\n")
  cat("\n")
} 