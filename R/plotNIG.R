#' Plot of LRR cummulative distribution function (CDF) of two models
#' 
#' This function plots the CDF of the LRR values and the expected CDF under Normal and Inverse Gaussian (NIG) distributions.
#' 
#' @param object An object of class MADloy from the \code{madloy} function.
#' @param tit Title of the plot.
#' 
#' @return A data table with the results for all samples in columns
#' @export
#' @examples
#' \dontrun{
#' plotNIG(resMADloy)}
plotNIG <- function(object, tit){
  if (!inherits(object, "MADloy"))
    stop("object must be of class 'madloy'")
  x <- MADloy:::getSummary(object)[,1]
  pars <- GeneralizedHyperbolic:::nigFit(x)
  pp <- pars$param
  pvals <- GeneralizedHyperbolic:::pnig(x, pp[1], pp[2], pp[3], pp[4])
  ss <- seq(min(x), max(x), len=100)
  plot(ecdf(x), xlab="LRR", main="")
  ss.norm <- pnorm(ss, mean=mean(x), sd=sd(x))
  ss.nig <- pnig(ss, pp[1], pp[2], pp[3], pp[4])
  ss.ecdf <- ecdf(x)(ss)
  lines(ss, ss.norm, col="red")
  lines(ss, ss.nig, col="blue")
  legend("topleft", c("Normal", "NIG"), col=c("red", "blue"), lty=1)
  title(tit)
  ans <- cbind(ss, ss.ecdf, ss.norm, ss.nig)
  colnames(ans) <- c("x", "ecdf", "norm", "nig")
  invisible(ans)
}
