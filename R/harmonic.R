loglikHarmonic <- function(x,par){ #par[1]=a, par[2]=m
  a <- par[1]
  m <- par[2]
  aux0 <- 1/(2*x*besselK(nu=0,x=a))
  aux1 <- exp(-(a/2)*(x/m + m/x))
  l <- sum(log(aux0) + log(aux1))
  if(is.infinite(l) || is.nan(l))
    l <- -1e6
  l
}

mleHarmonic <- function(x, a.ini=1, m.ini=1, ...){
  par.ini <- c(a.ini, m.ini)
  mle <- optim(par.ini, loglikHarmonic, x=x, method="L-BFGS-B", lower=c(1e-3, 1e-3), upper=c(Inf,Inf), control=list(fnscale=-1))
  mle
}

dHarmonic <- function(x, a, m){
  aux0 <- 1/(2*x*besselK(nu=0,x=a))
  aux1 <- exp(-(a/2)*(x/m + m/x))
  d <- aux0*aux1
  d
}

pHarmonic <- function(q, a, m){
  if(length(q)==1)
    p <- integrate(dHarmonic, lower=1e-5, upper=q, a=a, m=m)
  else{
    p <- NULL
    for(i in 1:length(q))
      p <- c(p, integrate(dHarmonic, lower=1e-5, upper=q[i], a=a, m=m)$value)
  }
  p
}

loglikNorm <- function(x, par){
  mu <- par[1]
  sd <- par[2]
  l <- sum(log(dnorm(x,mean=mu,sd=sd)))
  if(is.infinite(l) || is.nan(l))
    l <- -1e6
  l
}

