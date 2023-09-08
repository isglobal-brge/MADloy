#' Detection algorithm to detect Loss of Y events in MADloy data
#' 
#' @param x A MADloy object.
#' @param coef this determines how far the outliers are considered.
#' @param ... Other parameters.
#'   
#' @return An object of class 'LOY' that summarizes the LOY events detected in
#'   the analyzed samples
#' @export
#' @examples
#' \dontrun{
#' getLOY(resMADloy)}
getLOY <- function (x, coef=3, ...) 
{
  if (!inherits(x, "MADloy"))
    stop("x must be an object of class 'MADloy'")
  
## Not necessary (fixed JRG: 08/09/2023 - Maria's email)
#  xx <- getSummary(x)
#  target <- xx[,1]
#  reference <- xx[, 2]
  
  norm.lrr <- x$mLRRY$mLRRY

  findOutliers <- function(x, coef=coef) {
    upperq <- stats::quantile(x, 0.75, na.rm=TRUE)
    lowerq <- stats::quantile(x, 0.25, na.rm=TRUE)
    iqr <- stats::IQR(x, na.rm=TRUE)
    # we identify extreme outliers
    extreme.threshold.upper <- (iqr * coef) + upperq
    extreme.threshold.lower <- lowerq - (iqr * coef)
    ans <- x > extreme.threshold.upper | 
      x < extreme.threshold.lower
    ans
  }
  
  getThreshold <- function(x, ...)
  {
    den <- stats::density(x[!is.na(x)], bw="SJ", ...)
    m <- den$x[which(den$y==max(den$y))]
    xx <- x[x>=m]
    x.99 <- stats::quantile(xx, 0.99, na.rm=TRUE)
    ans <- m - x.99
    ans
  }
  
  outl <- findOutliers(norm.lrr, coef=coef)
  
  cl <- rep(NA, length(norm.lrr))
  cl[outl & norm.lrr>0] <- "XYY"
  cl[outl & norm.lrr<0] <- "LOY"
  cl[!outl] <- "normal"
  cl.f <- factor(cl, levels=c("normal", "LOY", "XYY"))
            
  tt <- getThreshold(norm.lrr)
  fosb <- stats::relevel(cut(norm.lrr, c(-Inf, tt, Inf),
              labels=c("LOY", "normal")),2)
  ans <- list()
  ans$res <- data.frame(MADloy = cl.f, Fosberg = fosb,
                    continous = norm.lrr)
  ans$data <- norm.lrr
  ans$par <- x$par
  attr(ans$data, "type") <- "LRR"
  class(ans) <- "LOY"
  
  ans
}
