#' @method plot MADseqLOY
#' @export
# TODO Set optimal png parameters to plot.
# TODO Optimize and improve performance.
plotIndNGS <-
  function(x, sample, label, file, readmax = 1000, logb = 10, rad = 1.5,
           deltatheta = 0, cexnum = 0.5, cexlab = 0.5, rado = 1.5, radi = 1, logged = TRUE,
           mlab = "Probes in Y region\n to detect LOY", filled = FALSE, ...) {
    if (!inherits(x, "MADseqLOY"))
      stop("object must be of class 'MADseqLOY' \n")
    if (missing(x))
      stop("A x object is required \n")
    if (missing(sample))
      stop("A name of the sample to be processed is required \n")
    
    ## Locate the index of the sample
    nn <- grep(sample, x$par$files)
    if (length(nn) == 0)
      stop("There are no samples in x object that matches '", sample, "'\n")
    ## Use the label if given. Else, use the name of the file.
    if (missing(label))
      label <- names(x[[1]][nn])
    
    ## Reading chrY for selected sample
    rawreads <- as.numeric(x$target[[nn]]$coverageAll[[1]])
    rawTargMedianCov <- x$target[[nn]]$medianTargetCoverage
    rawRef1MedianCov <- x$reference.1[[nn]]$medianTargetCoverage
    rawRef2MedianCov <- x$reference.2[[nn]]$medianTargetCoverage
    kbmax <- length(rawreads)
    lreads <-
      log(pmax(1, rawreads), logb)  ## Minimum value set to 1 to calculate log10
    ml <- min(lreads)
    lreads <- lreads - ml
    lreads <- lreads / (log(readmax, logb) - ml)
    reads <- rawreads / readmax
    
    pos <- 0:(kbmax - 1)
    thets <- pi / 2 + 2 * pi * pos / kbmax
    dxl <- (lreads + rad) * cos(thets)
    dyl <- (lreads + rad) * sin(thets)
    dx <- (reads + rad) * cos(thets)
    dy <- (reads + rad) * sin(thets)
    
    ## Target median coverage (chrY)
    TargMedianCov <- rawTargMedianCov / readmax  #Scaled
    
    lTargMedianCov <- log(rawTargMedianCov, logb)
    lTargMedianCov <- lTargMedianCov - ml
    lTargMedianCov <-
      lTargMedianCov / (log(readmax, logb) - ml)  #Log Scaled
    
    
    ## Reference 1 median coverage (chr21)
    Ref1MedianCov <- rawRef1MedianCov / readmax  #Scaled
    
    lRef1MedianCov <- log(rawRef1MedianCov, logb)
    lRef1MedianCov <- lRef1MedianCov - ml
    lRef1MedianCov <-
      lRef1MedianCov / (log(readmax, logb) - ml)  #Log Scaled
    
    ## Reference 2 median coverage (chr22)
    
    Ref2MedianCov <- rawRef2MedianCov / readmax  #Scaled
    
    lRef2MedianCov <- log(rawRef2MedianCov, logb)
    lRef2MedianCov <- lRef2MedianCov - ml
    lRef2MedianCov <-
      lRef2MedianCov / (log(readmax, logb) - ml)  #Log Scaled
    
    ## Check if png
    if (!missing(file))
      grDevices::png(
        filename = file, width = 150, height = 150, units = "mm", res = 300,
        pointsize = 20 
      )
    
    ## Shift text by small angle to make it centred on the region it's labelling...
    offset <- 2 * pi * deltatheta / 360
    op <- graphics::par(mai = c(0, 0, 0, 0), bg = "white")
    
    ## Draw a circle
    graphics::plot(
      NULL, xlim = c(-2.7, 2.7), ylim = c(-2.7, 2.7), axes = FALSE, ann = FALSE
    )
    
    
    ## Trace scale
    rads <- (10 ^ seq(0, 3))
    lrads <-
      (log(rads, logb) - log(1, logb)) / (log(readmax, logb) - log(1, logb))
    radslin <- seq(0, readmax, 500)
    if (!logged) {
      rads <- radslin
      lrads <- radslin / readmax
    }
    for (j in 1:length(rads)) {
      rd <- lrads[j]
      ecolitk::linesCircle(rado + rd, lty = 2)
      graphics::text(0, 1.025 * rado + rd, rads[j], cex = 0.5)
    }
    
    ## Draw histogram
    tlist <- clist <- ltylist <- NULL
    tlist <- c(tlist, "Targets", "Target region coverage")
    clist <- c(clist, "black", "blue")
    ltylist <- c(ltylist, 1)
    
    if (filled) {
      if (logged) {
        graphics::polygon(dxl, dyl, col = "blue", border = NA)
      } else {
        graphics::polygon(dx, dy, border = NA, col = "blue")
      }
    } else {
      if (logged) {
        graphics::points(dxl, dyl, type = "l", col = "blue", lwd = 2)
      } else {
        graphics::points(dx, dy, type = "l", col = "blue", lwd = 2)
      }
    }
    plotrix::draw.circle(
      0, 0, rado, nv = 100, border = "white", col = "white", lty = 1, lwd = 1
    )
    
    
    ## Draw probes locations
    
    
    probeRanges <-
      IRanges::ranges(x$par$targets[IRanges::overlapsAny(x$par$targets, x$par$target.region)])
    
    for (i in 1:length(probeRanges)) {
      end = end(probeRanges)[i]
      start = start(probeRanges)[i]
      ecolitk::polygonArc(
        theta0 = pi / 2 + 2 * pi * start / kbmax, theta1 = pi / 2 + 2 * pi *
          end / kbmax, radius.in = 1.25, radius.out = rado, edges = 1000, col = "black",
        border = "black", lwd = 1
      )
    }
    
    ## Draw Medians
    ecolitk::linesCircle(rado + lTargMedianCov, lty = 2, lwd = 3, col = "green")
    tlist <-
      c(tlist, paste0("Target median coverage= ", rawTargMedianCov))
    clist <- c(clist, "green")
    ltylist <- c(ltylist, 2)
    
    ecolitk::linesCircle(rado + lRef1MedianCov, lty = 2, lwd = 3, col = "red")
    tlist <-
      c(tlist, paste0("Reference 1 median coverage = ", rawRef1MedianCov))
    clist <- c(clist, "red")
    ltylist <- c(ltylist, 2)
    
    ecolitk::linesCircle(rado + lRef2MedianCov, lty = 2, lwd = 3, col = "yellow")
    tlist <-
      c(tlist, paste0("Reference 2 median coverage = ", rawRef2MedianCov))
    clist <- c(clist, "yellow")
    ltylist <- c(ltylist, 2)
    
    ## Draw coordinate numbers
    kbs <- seq(0, kbmax, 1e+06)
    trad <- 0.75 * radi
    for (kb in kbs) {
      theta <- pi / 2 + 2 * pi * (kb / max(pos))
      theta <- theta + offset
      ttext <- pi + theta
      strw <- graphics::strwidth(paste(kb / 1000, " mb", sep = ""), cex = cexnum)
      dist <- trad + strw
      graphics::text(
        dist * cos(theta), dist * sin(theta), paste(kb / 1e+06, " mb", sep = ""),
        col = "black", cex = 1 * cexnum, srt = 360 * ttext / (2 * pi), adj = c(0,
                                                                               0)
      )
    }
    
    ## Add Text in the plot
    graphics::text(0, 0, mlab, cex = 1)
    
    ## Add Title
    graphics::text(0, 2.75, label, cex = 1.05)
    
    ## Add Legend
    graphics::legend(
      "bottomright", legend = tlist, lty = ltylist, lwd = 2, col = clist, bty = "n",
      cex = 0.75
    )
    
    ## Parameters
    graphics::par(op)
    
    ## Close png
    if (!missing(file))
      grDevices::dev.off()
  }
