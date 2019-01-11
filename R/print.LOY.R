#' @export
print.LOY <- function(x, method="MADloy", ...) {
  
  if (method=="MADloy")
    data <- x$MADloy
  else if (method=="Fosberg")
    data <- x$Fosberg
  else
    stop("method should be either 'MADloy' or 'Fosberg'")
       
  cat("\n")
  cat("Object of class LOY \n")
  cat("---------------------------- \n")
  cat("Number of normal samples:", sum(data == "normal", 
                                       na.rm=TRUE), "\n")
  cat("Number of LOY:", sum(data == "LOY", 
                            na.rm=TRUE), "\n")
  cat("Number of XYY:", sum(data == "XYY", 
                            na.rm=TRUE), "\n")
  cat("Number of samples do not pass QC:", sum(is.na(data)))
  cat("\n")
}
