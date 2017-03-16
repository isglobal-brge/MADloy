plotNIG <- function(x, tit){
  pars <- nigFit(x)
  pp <- pars$param
  pvals <- pnig(x, pp[1], pp[2], pp[3], pp[4])
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
