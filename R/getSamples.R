#' @export
getSamples <- function(x, type="LOY", ...) {
  if(!inherits(x, c("MADloyBdev", "MADloy", "LOY")))
    stop("x should be an object of class 'MADloyBdev', 'MADloy' or 'LOY'")
  class <- as.factor(x$class)
  names(class) <- tools:::file_path_sans_ext(x$par$files)
  ans <- names(class)[class%in%type]
  ans
}
