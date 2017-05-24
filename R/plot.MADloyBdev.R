#' @method plot MADloyBdev
#' @export
plot.MADloyBdev <- function(x, ...) {
  sel <- x$class$adjusted_p != "balancedpq"
  filled <- rep(1, sum(sel))
  filled[(x$class$adjusted_p != "balancedpq")[sel]] <- 20
  
  gplots::plotCI(x = unlist(data.frame(t(sapply( x$Bdev, "[[", "p")))$Pl[sel]), 
         uiw = unlist(data.frame(t(sapply( x$Bdev, "[[", "p")))$Plsd[sel]), 
         lty = 2, 
         xaxt ="n", 
         pch=filled, 
         gap = 0, 
         ylim=c(0, 3), 
         minbar = 0, 
         maxbar = 3, 
         main="Ploidy trimmed mean and deviation in p and q arms", 
         ylab="Ploidy", 
         xlab="getLOY classification")
  
  par(new=TRUE)
  
  gplots::plotCI(x = unlist(data.frame(t(sapply( x$Bdev, "[[", "q")))$Pl[sel]), 
         uiw = unlist(data.frame(t(sapply( x$Bdev, "[[", "q")))$Plsd[sel]), 
         lty = 2, 
         xaxt ="n", 
         pch=filled, 
         gap = 0, 
         ylim=c(0, 3), 
         minbar = 0, 
         maxbar = 3, 
         col="red", 
         ylab="", 
         xlab="")
  
  axis(1, 
       at=c(1:sum(sel)), 
       labels=FALSE)
  
  text(x=c(1:sum(sel)), 
       y=par()$usr[3]-0.125,
       labels=x$class$orig[sel], 
       srt=45, 
       adj=1, 
       xpd=TRUE, 
       cex=0.75)
  
  legend("topright", 
         legend = c("p arm", "q arm"), 
         pch = c(20, 20), 
         col=c("black", "red"))
} 
