#' Detection algorithm to detect Loss of Y events in MADloy or MADseqLOY data
#' 
#' @param object A MADloy or MADseqLOY object.
#' @param ref One of the reference chromosomes used in the \code{MADloy} or \code{MADseqLOY}
#'   functions.
#' @param k Number of expected classification groups
#' @param ... Other parameters
#'   
#' @return An object of "LOY" class that summarizes the LOY events detected in
#'   the analyzed samples
#' @export
#' @examples
#' \dontrun{
#' getLOY(resMADseqLOY)
#' getLOY(resMADloy)}
getLOY <- function(object, ref = "22", k , cutoff= 0.90, arbvar=FALSE, pval.harm=0.05, ...) {
  
  if (inherits(object, "MADseqLOY") | inherits(object, "MADloy")) {
    x <- MADloy:::getMedian(object)
  } else {
    x <- object
  }
  target <- x[, 1]
  posRef <- grep(ref, colnames(x))
  if (length(posRef) == 0) 
    stop(paste0("Chromosome reference ", ref, " is not valid. Check arguments 'ref.region.1' and 'ref.region.2' in MADseqLOY funcion"))
  reference <- x[, posRef]
  
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
      stop("Model does not converge: change initial parameters or run again the function")

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
      stop("Model does not converge: change initial parameters or run again the function")
     p.test <- 1
    }

    cl <- apply(ans$posterior, 1, which.max)

    if (p.test >= 0.01) {
    ratio <- xx[, 1]/xx[, 2]
    tt <- aggregate(ratio ~ as.factor(cl), FUN=mean)
    o <- order(tt[,2])
     alt <- which.max(abs(tt[,2] - 1))
     if (tt[alt,2] > 1) {
       labs <- c("normal", "GAIN") }
     else {
       labs <- c("LOY", "normal") }
    cl <- factor(cl, labels = labs[o])     
    }
    else {
     ratio <- xx[, 1]/xx[, 2]
     tt <- aggregate(ratio ~ as.factor(cl), FUN=mean)
     o <- order(tt[,2])
     cl <- factor(cl, labels = c("LOY", "normal", "GAIN")[o])     
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
  } else {
    xx <- cbind(target, reference)
    y <- target-reference
    ref <- reference
    pars <- nigFit(ref, trace=FALSE)
    pp <- pars@fit$estimate
    pvals <- pnig(y, pp[1], pp[2], pp[3], pp[4])
    
    cl <- ifelse(pvals > pval.harm, "normal", "altered")
    cl[cl=="altered" & log(y) > 0] <- "GAIN"
    cl[cl=="altered" & log(y) < 0] <- "LOY"
    ans <- list(class = cl, prob = pvals, data = xx)
    attr(ans, "type") <- "LRR"
    class(ans) <- "LOY"
  }
  ans
} 
