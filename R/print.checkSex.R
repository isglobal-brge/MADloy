#' @export
print.checkSex <- function(x, ...) {
    cat("Object of class checkSex \n")
    cat("---------------------------- \n")
    cat("Number of processed samples:", length(x$class), "\n")
    cat("Number of female classified samples: ", sum(x$class == "FEMALE"), "\n")
    cat("Number of male classified samples: ", sum(x$class == "MALE"), "\n")
    cat("\n")
}
