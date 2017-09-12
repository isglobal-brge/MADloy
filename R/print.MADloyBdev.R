#' @export
print.MADloyBdev <- function(x, ...) {
    cat("Object of class MADloyBdev \n")
    cat("---------------------------- \n")
    print(table(LRRClassification=as.character(x$class$orig), BdevClassification=x$class$class))
    cat("\n")
}
