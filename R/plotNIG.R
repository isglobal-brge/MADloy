plotNIG <- function(x, tit){
  pars <- nigFit(x)
  pp <- pars$param
  pvals <- pnig(x, pp[1], pp[2], pp[3], pp[4])
  ss <- seq(min(x), max(x), len=100)
  plot(ecdf(x), main="")
  lines(ss, pnorm(ss, mean=mean(x), sd=sd(x)), col="red")
  lines(ss, pnig(ss, pp[1], pp[2], pp[3], pp[4]), col="blue")
  legend("topleft", c("Normal", "NIG"), col=c("red", "blue"), lty=1)
  title(tit)
}
