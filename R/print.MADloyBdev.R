#' @export
print.MADloyBdev <- function(x, ...) {
    cat("Object of class MADloyBdev \n")
    cat("---------------------------- \n")
    class <- as.factor(x$class)
    names(class) <- tools:::file_path_sans_ext(x$par$files)
    print(table(class))
    cat("\n")
}
