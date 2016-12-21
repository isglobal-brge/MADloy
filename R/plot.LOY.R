#' @method plot LOY
#' @export
plot.LOY <- function(x, labels, colors = c("red", "blue", "darkgreen"), pos.leg = "bottomleft", ...) {
  data <- x$data
  if (missing(labels)) {
    labels <- rownames(data)
  }
  
  xx <- data[, 1]
  yy <- data[, 2]
  
  nclass <- length(unique(x$class))

  if (nclass==2) {
   loy <- x$class == "LOY"
   gain <- x$class == "GAIN"
   if (sum(loy) > sum(gain)) { 
     mycol <- ifelse(loy, colors[1], colors[2])
     leg.lab <- c("LOY", "normal") 
     col.lab <- colors[1:2]
     alt <- loy }
   else {
     mycol <- ifelse(gain, colors[3], colors[2])
     leg.lab <- c("normal", "GAIN")
     col.lab <- colors[2:3]
     alt <- gain}  
  
  }
  if (nclass==3) {
   loy <- x$class == "LOY"
   gain <- x$class == "GAIN"
   mycol <- ifelse(loy, colors[1], ifelse(gain, colors[3],colors[2]))
   leg.lab <- c("LOY", "normal", "GAIN")
   col.lab <- colors
   alt <- x$class%in%c("LOY", "GAIN")
  }  

  d <- xx - yy
  ss <- 1:nrow(data)
  plot.default(ss, d, type = "n", xlab = "Individuals", ylab = "Mean coverage difference (Reference - mY)")
  points(ss, d, col = mycol, pch = 16, ...)
  legend(pos.leg, leg.lab, pch = 16, col = col.lab)
  text(ss[alt], jitter(d[alt]), labels[alt], cex = 0.8, adj = 0)
 # if (nclass==3)
 #  text(ss[gain], jitter(d[gain]), labels[gain], cex = 0.8, adj = 0)
} 




plot.LOY2 <- function(x, labels, colors = c("red", "blue", "darkgreen"), pos.leg = "bottomleft", ...) {
  data <- x$data
  if (missing(labels)) {
    labels <- rownames(data)
  }
  
  xx <- data[, 1]
  yy <- data[, 2]
  
  nclass <- length(unique(x$class))

  if (nclass==2) {
   loy <- x$class == "LOY"
   gain <- x$class == "GAIN"
   if (sum(loy) > sum(gain)) { 
     mycol <- ifelse(loy, colors[1], colors[2])
     leg.lab <- c("LOY", "normal") 
     col.lab <- colors[1:2]
     alt <- loy }
   else {
     mycol <- ifelse(gain, colors[3], colors[2])
     leg.lab <- c("normal", "GAIN")
     col.lab <- colors[2:3]
     alt <- gain}  
  
  }
  if (nclass==3) {
   loy <- x$class == "LOY"
   gain <- x$class == "GAIN"
   mycol <- ifelse(loy, colors[1], ifelse(gain, colors[3],colors[2]))
   leg.lab <- c("LOY", "normal", "GAIN")
   col.lab <- colors
   alt <- x$class%in%c("LOY", "GAIN")
  }  

  d <- xx - yy
  ss <- 1:nrow(data)
  plot.default(yy, xx, type = "n", ylab = "Median coverage (male-Y region)", xlab = "Mean coverage (Reference)")
  points(yy, xx, col = mycol, pch = 16, ...)
  legend(pos.leg, leg.lab, pch = 16, col = col.lab)
  text(yy[alt], jitter(xx[alt]), labels[alt], cex = 0.8, adj = 0)
  if (nclass==3)
   text(yy[gain], jitter(xx[gain]), labels[gain], cex = 0.8, adj = 0)
} 
