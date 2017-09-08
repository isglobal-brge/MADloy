#' Detection algorithm to detect Loss of Y events in MADloy or MADseqLOY data
#' 
#' @param object A MADloy or MADseqLOY object.
#' @param k Number of groups. Only necessary in NGS data.
#' @param cutoff cutoff value. Only necessari in NGS data.
#' @param pval.sig pval.sig p-value treshold to be used in the classification test.
#' @param ... Other parameters.
#'   
#' @return An object of class 'LOY' that summarizes the LOY events detected in
#'   the analyzed samples
#' @export
#' @examples
#' \dontrun{
#' getLOY(resMADseqLOY)
#' getLOY(resMADloy)}
getLOY <- function(object, pval.sig, k, cutoff, ...) {
    
    if (inherits(object, "MADseqLOY") | inherits(object, "MADloy")) {
        x <- getSummary(object)
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
            while (control < 10 & (inherits(ans2, "try-error") | inherits(ans3, "try-error"))) {
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
            p.test <- stats::pchisq(-2 * (ans2$loglik - ans3$loglik), df.diff, lower = FALSE)
            if (p.test < 0.01) {
                ans <- ans3
            } else {
                ans <- ans2
            }
        } else {
            control <- 0
            ans <- list()
            class(ans) <- "try-error"
            while (control < 10 & inherits(ans, "try-error")) {
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
            tt <- stats::aggregate(ratio ~ as.factor(cl), FUN = mean)
            o <- order(tt[, 2])
            alt <- which.max(abs(tt[, 2] - 1))
            if (tt[alt, 2] > 1) {
                labs <- c("normal", "XYY")
            } else {
                labs <- c("LOY", "normal")
            }
            cl <- factor(cl, labels = labs[o])
        } else {
            ratio <- xx[, 1]/xx[, 2]
            tt <- stats::aggregate(ratio ~ as.factor(cl), FUN = mean)
            o <- order(tt[, 2])
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
    } else {
        
        if (missing(pval.sig)) 
            pval.sig <- 0.05/length(target)
        
        norm.lrr <- object$mLRRY
        reference <- unlist(lapply(egcut$reference, "[[", "summary"))
        reference.qc <- reference[!is.na(norm.lrr)]
        
        sds <- unlist(lapply(egcut$reference, "[[", "sd"))
        
        
        pars <- GeneralizedHyperbolic::nigFit(reference.qc)
        # pars <- fBasics:::nigFit(ref, trace=FALSE) pp <- pars@fit$estimate
        pp <- pars$param
        
        ff <- function(x, param) {
            if (x > 0 & !is.na(x)) 
                ans <- try(GeneralizedHyperbolic::pnig(x, param[1], param[2], param[3], 
                  param[4], lower = FALSE), TRUE) else ans <- try(GeneralizedHyperbolic::pnig(x, param[1], param[2], param[3], 
                param[4], lower = TRUE), TRUE)
            if (is.na(x))
                ans <- NA
            if (inherits(ans, "try-error")) 
                ans <- NA
            return(ans)
        }
        
        pvals <- sapply(norm.lrr, ff, param = pp)
        
        threshold <- stats::median(sds[!sds > 2 * mean(sds)])
        cl <- ifelse(pvals > pval.sig | abs(norm.lrr) < threshold, "normal", "altered")
        cl[cl == "altered" & norm.lrr > 0] <- "XYY"
        cl[cl == "altered" & norm.lrr < 0] <- "LOY"
        cl.o <- factor(ans$class, levels=c("normal", "LOY", "XYY"))
        par <- object$par
        par$offset <- offset
        par$pval.sig <- pval.sig
        ans <- list(class = cl.o, prob = pvals, data = norm.lrr, ref = reference, par = par)
        attr(ans, "type") <- "LRR"
        class(ans) <- "LOY"
    }
    ans
}
