processBam <- function(file, targets, subset) {
    targets <- targets[IRanges::overlapsAny(targets, subset)]
    bamConn <- Rsamtools::BamFile(file)
    GenomeInfoDb::seqlevelsStyle(subset) <- GenomeInfoDb::seqlevelsStyle(bamConn)
    if (BiocGenerics::end(subset) == 3e+08) 
        BiocGenerics::end(subset) <- GenomeInfoDb::seqlengths(bamConn)[GenomeInfoDb::seqlevels(subset)]
    param <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE), 
        what = c("qname", "pos", "qwidth", "rname"), which = subset)
    open(bamConn)
    aln <- Rsamtools::scanBam(file = bamConn, param = param)[[1]]
    reads <- with(aln, GenomicRanges::GRanges(IRanges::IRanges(start = pos, width = qwidth), 
        ID = qname, seqnames = rname))
    close(bamConn)
    GenomeInfoDb::seqlevelsStyle(reads) <- "UCSC"
    GenomeInfoDb::seqlevelsStyle(subset) <- "UCSC"
    reads <- methods::as(reads, "RangedData")[as.character(GenomeInfoDb::seqnames(subset))]
    targets <- methods::as(targets, "RangedData")[as.character(GenomeInfoDb::seqnames(subset))]
    coverage <- TEQC::coverage.target(reads = reads, targets = targets)
    coverage$coverageAll <- NULL
    coverage$medianTargetCoverage <- stats::median(coverage$coverageTarget[[1]])
    return(coverage)
}
