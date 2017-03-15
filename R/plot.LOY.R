#' @method plot LOY
#' @export
plot.LOY <- function(x, labels, colors = c("red", "blue", "darkgreen"), pos.leg = "bottomleft", ...) {
  data <- x$data
  if (missing(labels)) {
    if(attr(x, "type")=="Coverage")
     labels <- rownames(data)
    if(attr(x, "type")=="LRR")
      labels <- names(data)
  }
  
  
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

  
  if(attr(x, "type")=="Coverage") {
   ss <- 1:nrow(data)
   xx <- data[, 1]
   yy <- data[, 2]  
   d <- xx - yy
   plot.default(ss, d, type = "n", xlab = "Individuals", ylab = "Mean coverage difference (mY region - Reference)")
   points(ss, d, col = mycol, pch = 16, ...)
  }
  if(attr(x, "type")=="LRR") {
    ss <- 1:length(data)
    d <- x$data
    plot.default(ss, d, type = "n", xlab = "Individuals", 
                 ylab = "Mean normalized LRR in mY region")
    points(ss, d, col = mycol, pch = 16, ...)
  }
  legend(pos.leg, leg.lab, pch = 16, col = col.lab)
  alt <- x$class%in%c("LOY", "GAIN")
  wordcloud::textplot(x = ss[alt], y = d[alt], words = tools::file_path_sans_ext(labels[alt]), 
                      cex = 0.8, new = FALSE, xlim=c(min(ss[alt]), max(ss[alt])), ylim=c(min(d[alt]), max(d[alt])))
 # if (nclass==3)
 #  text(ss[gain], jitter(d[gain]), labels[gain], cex = 0.8, adj = 0)
} 

