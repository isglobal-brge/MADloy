#' Get probability matrix of three status (LOY, normal, XYY)
#' 
#' @param x an object of class 'LOY'
#' @return a matrix of probabilities
#' 
#' @export
#' 
#' 
getProb <- function(x) {
  n <- length(x$class)
  idx <- sapply(x$class, function(x) which(x== c("LOY", "normal", "XYY")))
  ans <- matrix(0, ncol=3, nrow=n)
  
  for (i in 1:n) {
    if (idx[i]==1)
      ans[i, 1:2] <- c(1-x$prob[i], x$prob[i]) 
    if (idx[i]==2 & x$data[i]>=0) 
      ans[i, 2:3] <- c(1-x$prob[i], x$prob[i]) 
    if (idx[i]==2 & x$data[i]<0) 
      ans[i, 1:2] <- c(x$prob[i], 1-x$prob[i])
    if (idx[i]==3)
      ans[i, 2:3] <- c(x$prob[i], 1-x$prob[i]) 
  }
  ans
}
