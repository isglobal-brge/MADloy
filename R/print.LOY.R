#' @export
print.LOY <- function(x, ...) {
    p <- unlist(x$estim[1])
    # se <- x$se[1] ci <- c(p - 1.96 * se, p + 1.96 * se)
    
    cat("\n")
    cat("Object of class LOY \n")
    cat("---------------------------- \n")
    cat("Number of normal samples:", sum(x$class == "normal"), "\n")
    cat("Number of LOY:", sum(x$class == "LOY"), "\n")
    cat("Number of XYY:", sum(x$class == "XYY"), "\n")
    cat("\n")
    
    # TODO adapt to MADloy processed objects cat('Estimated proportion of LOY (CI
    # 95%):', round(p * 100, 1), '(', round(ci[1] * 100, 1), '-', round(ci[2] * 100,
    # 1), ') \n')
}
