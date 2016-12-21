read.targets <- function(targetsfile, chrcol = 1, startcol = 2, endcol = 3, zerobased = TRUE, 
  sep = "\t", skip = 1, header = FALSE, quiet, ...) {
  
  n <- ncol(data.table::fread(targetsfile, sep = sep, skip = 60, nrows = 1))
  colclasses <- rep("NULL", n)
  colclasses[chrcol] <- "character"
  colclasses[c(startcol, endcol)] <- "integer"
  new_colnames <- rep("NULL", 3)
  new_colnames[c(chrcol, startcol, endcol)] <- c("chrcol", "startcol", "endcol")
  dat <- data.table::fread(targetsfile, sep = sep, colClasses = colclasses, skip = skip, 
    header = header, ...)
  setnames(dat, colnames(dat), new_colnames)
  ir <- IRanges::IRanges(start = dat$startcol, end = dat$endcol)
  if (zerobased) 
    BiocGenerics::start(ir) <- BiocGenerics::start(ir) + 1
  gr <- GenomicRanges::GRanges(seqnames = dat$chrcol, ranges = ir)
  n.targ <- length(gr)
  gr <- GenomicRanges::reduce(gr)
  n.targ.no <- length(gr)
  if (!quiet) {
    if (n.targ == n.targ.no) {
      message(paste("Read", n.targ, "(non-overlapping) exome target regions\n"))
    } else {
      message(paste("Read", n.targ, "exome target regions in total, which are collapsed to", 
        n.targ.no, "non-overlapping exome target regions\n"))
    }
  }
  return(gr)
} 
