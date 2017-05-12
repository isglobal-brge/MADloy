#' Detection algorithm to detect Loss of Y events in MADloy or MADseqLOY data
#' 
#' @param object A MADloy or MADseqLOY object.
#' @param offset Offset value of SNP array data in the LRR of chromosome Y. That is, value to
#' guarantee that mean LRR at chrmosome Y is centered at 0.  Default 0 since LRR at 
#' m-LRR region is expected to be centered at O. 
#' @param k Number of groups. Only necessary in NGS data.  
#' @param ... Other parameters.
#'   
#' @return An object of class "LOY" that summarizes the LOY events detected in
#'   the analyzed samples
#' @export
#' @examples
#' \dontrun{
#' getLOY(resMADseqLOY)
#' getLOY(resMADloy)}
getLOY <- function(object, offset=0, pval.sig=0.05, ...) {
  
  if (inherits(object, "MADseqLOY") | inherits(object, "MADloy")) {
    x <- MADloy:::getSummary(object)
  } else {
    x <- object
  }
  target <- x[, 1]
  reference <- x[, 2]
  
  if (inherits(object, "MADseqLOY")) {
    xx <- cbind(target, reference)
    if (missing(k)) {
     control <- 0
     ans2 <- ans3 <- list()
     class(ans2) <- class(ans3) <- "try-error"
     while (control<10 & (inherits(ans2, "try-error") | inherits(ans3, "try-error"))) {
       if (inherits(ans2, "try-error"))
        ans2 <- try(mixtools::mvnormalmixEM(xx, k = 2, ...), TRUE)
       if (inherits(ans3, "try-error"))
       ans3 <- try(mixtools::mvnormalmixEM(xx, k = 3, ...), TRUE) 
       control <- control + 1
     } 
     if (inherits(ans2, "try-error") | inherits(ans3, "try-error"))
      stop("Model does not converge: change initial parameters or run aXYY the function")

     df2 <- length(unique(unlist(ans2$estim)))
     df3 <- length(unique(unlist(ans3$estim)))
     df.diff <- df3 - df2
     p.test <- pchisq(-2*(ans2$loglik - ans3$loglik), df.diff, lower=FALSE)
     if (p.test < 0.01) {
       ans <- ans3 }
     else {
       ans <- ans2
     }
    }
    else {
     control <- 0
     ans <- list()
     class(ans) <- "try-error"
     while(control < 10 & inherits(ans, "try-error")) {
      ans <- try(mixtools::mvnormalmixEM(xx, k = k, ...), TRUE)
      control <- control + 1 
     }
     if (inherits(ans, "try-error"))
      stop("Model does not converge: change initial parameters or run aXYY the function")
     p.test <- 1
    }

    cl <- apply(ans$posterior, 1, which.max)

    if (p.test >= 0.01) {
    ratio <- xx[, 1]/xx[, 2]
    tt <- aggregate(ratio ~ as.factor(cl), FUN=mean)
    o <- order(tt[,2])
     alt <- which.max(abs(tt[,2] - 1))
     if (tt[alt,2] > 1) {
       labs <- c("normal", "XYY") }
     else {
       labs <- c("LOY", "normal") }
    cl <- factor(cl, labels = labs[o])     
    }
    else {
     ratio <- xx[, 1]/xx[, 2]
     tt <- aggregate(ratio ~ as.factor(cl), FUN=mean)
     o <- order(tt[,2])
     cl <- factor(cl, labels = c("LOY", "normal", "XYY")[o])     
    }
    
    labs <- rownames(xx)
    probs <- ans$posterior
    rownames(probs) <- labs

    mprob <- apply(probs, 1, max)

    cl[cl != "normal" & mprob < cutoff] <- "normal"

    estim <- list(pi = ans$lambda, mu = ans$mu, sigma = ans$sigma)
     
   
    ans <- list(class = cl, estim = estim, prob = probs, logLik = ans$loglik, 
      data = xx)
    attr(ans, "type") <- "Coverage"
    class(ans) <- "LOY"
  } 
  else {
    xx <- cbind(target, reference)
    norm.lrr <- target - reference + offset
    pars <- GeneralizedHyperbolic:::nigFit(reference)
#    pars <- fBasics:::nigFit(ref, trace=FALSE)
#    pp <- pars@fit$estimate
    pp <- pars$param
    
    ff <- function(x, param){
      if (x>=0.1)
        ans <- GeneralizedHyperbolic:::pnig(x, param[1], param[2], 
                                            param[3], param[4], lower=FALSE)
      else
        ans <- GeneralizedHyperbolic:::pnig(x, param[1], param[2],
                                                param[3], param[4], lower=TRUE)
    }
    pvals <- sapply(norm.lrr, ff, param=pp)
    
    cl <- ifelse(pvals > pval.sig | abs(norm.lrr)<0.15, "normal", "altered")
    cl[cl=="altered" & norm.lrr > 0] <- "XYY"
    cl[cl=="altered" & norm.lrr < 0] <- "LOY"
    par <- object$par
    par$offset <- offset
    par$pval.sig <- pval.sig
    ans <- list(class = cl, prob = pvals, data = norm.lrr, par = par)
    attr(ans, "type") <- "LRR"
    class(ans) <- "LOY"
  }
  ans
} 
