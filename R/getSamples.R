#' Return the samples with the chosen classification
#' 
#' @param x A MADloy, LOY or MADloyBdev object.
#' @param type A string with the chosen classification.
#' @return A vector of the samples with the chosen classification.
#' @export
#' @examples
#' \dontrun{
#' getSamples(resMADloy)}

getSamples <- function(x, type="LOY") {
  if(!inherits(x, c("MADloyBdev", "MADloy", "LOY")))
    stop("x should be an object of class 'MADloyBdev', 'MADloy' or 'LOY'")
  class <- as.factor(x$calling$MADloy)
  names(class) <- tools::file_path_sans_ext(x$par$files)
  ans <- names(class)[class%in%type]
  ans
}
