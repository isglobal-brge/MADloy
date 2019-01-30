#' Detection algorithm to detect Loss of Y events in MADloy or MADseqLOY data
#' 
#' @param object A MADloy or MADseqLOY object.
#' @param coef this determines how far the outliers are considered.
#' @param ... Other parameters.
#'   
#' @return An object of class 'LOY' that summarizes the LOY events detected in
#'   the analyzed samples
#' @export
#' @examples
#' \dontrun{
#' getLOY(resMADseqLOY)
#' getLOY(resMADloy)}
getLOY <- function (object, coef=3, ...) 
{
  if (inherits(object, "MADseqLOY") | inherits(object, "MADloy")) {
    if (inherits(object, "MADseqLOY"))
        stop("method should be implemented (bi-normal or outlier)")
    x <- MADloy:::getSummary(object)
  }
  else {
    x <- object
  }
  target <- x[, 1]
  reference <- x[, 2]
  
  norm.lrr <- object$mLRRY

  findOutliers <- function(x, coef=coef) {
    upperq <- quantile(x, 0.75, na.rm=TRUE)
    lowerq <- quantile(x, 0.25, na.rm=TRUE)
    iqr <- IQR(x, na.rm=TRUE)
    # we identify extreme outliers
    extreme.threshold.upper <- (iqr * coef) + upperq
    extreme.threshold.lower <- lowerq - (iqr * coef)
    ans <- x > extreme.threshold.upper | 
      x < extreme.threshold.lower
    ans
  }
  
  getThreshold <- function(x, ...)
  {
    den <- density(x[!is.na(x)], bw="SJ", ...)
    m <- den$x[which(den$y==max(den$y))]
    xx <- x[x>=m]
    x.99 <- quantile(xx, 0.99, na.rm=TRUE)
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
  fosb <- relevel(cut(norm.lrr, c(-Inf, tt, Inf),
              labels=c("LOY", "normal")),2)
  ans <- list()
  ans$res <- data.frame(MADloy = cl.f, Fosberg = fosb,
                    continous = norm.lrr)
  ans$data <- object$mLRRY
  ans$par <- object$par
  attr(ans$data, "type") <- "LRR"
  class(ans) <- "LOY"
  
  ans
}
