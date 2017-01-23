#' Check the loss of Y events in a given folder of .bam files
#' 
#' MADseqLOY check all the bams in a given folder, by using the target regions
#' of the technology used to sequence the sample. The coverage is computed for
#' the region chrY:2694521-59034049 (hg19/GRCh37) and another given by the query argument.
#' 
#' @param files A folder path where the .bam files are or a vector of file
#'   paths. This function searches in recursive folders.
#' @param exomeTargets A file path to a .bed file with the regions targeted by
#'   the sequencing method used.
#' @param target.region The chromosome or region to be compared with the other
#'   regions. By default is the region chrY:2694521-59034049 but it can be
#'   changed.
#' @param ref.region.1 First chromosome or region to be compared with the Y
#'   region in UCSC style (i.e. "chr21" or "chr21:1000-10000").
#' @param ref.region.2 Second chromosome or region to be compared with the Y
#'   region in UCSC style (i.e. "chr22" or "chr22:1000-10000").
#' @param mc.cores umber of cores to be used with this function. If there are
#'   more cores than samples, the number of cores will be limited to the number
#'   of samples. By default set to 1.
#' @param quiet Should the function not inform about the status of the process.
#'   By default is FALSE.
#' @param skip Number of lines to be skipped in the targets file if necessary.
#'   By default is set to 1.
#' @param ... Other parameters of the read.targets function.
#'   
#' @return The returned value is an object of "MADseqLOY" class with four lists,
#'   one regarding the target region coverage called "target", two regarding the
#'   reference regions to be compared against called "reference.1" and
#'   "reference.2" and finally one with the parameters used in the query.
#' @export
#' @examples
#' \dontrun{
#' madseqloy(files=bamFile, reference=targetFile, skip=0)}


madseqloy <- function (files, exomeTargets, target.region = "chrY:2694521-59034049", 
                       ref.region.1 = "chr21", ref.region.2 = "chr22", mc.cores, 
                       quiet = FALSE, skip = 2, ...) {
  if (missing(target.region)) 
    message("Targeted region set to chrY:2694521-59034049 by default\n")
  queryA <- unlist(strsplit(x = target.region, split = "[:, -]", 
                            perl = T))
  if (is.na(queryA[2]) | is.na(queryA[3])) {
    subsetA <- GenomicRanges::GRanges(seqnames = queryA[1], 
                                      ranges = IRanges::IRanges(start = 1, end = 3e+08))
  }
  else {
    subsetA <- GenomicRanges::GRanges(seqnames = queryA[1], 
                                      ranges = IRanges::IRanges(start = as.numeric(queryA[2]), 
                                                                end = as.numeric(queryA[3])))
  }
  if (missing(ref.region.1)) 
    message("Reference region 1 set to chr21 by default\n")
  queryB <- unlist(strsplit(x = ref.region.1, split = "[:, -]", 
                            perl = T))
  if (is.na(queryB[2]) | is.na(queryB[3])) {
    subsetB <- GenomicRanges::GRanges(seqnames = queryB[1], 
                                      ranges = IRanges::IRanges(start = 1, end = 3e+08))
  }
  else {
    subsetB <- GenomicRanges::GRanges(seqnames = queryB[1], 
                                      ranges = IRanges::IRanges(start = as.numeric(queryB[2]), 
                                                                end = as.numeric(queryB[3])))
  }
  if (missing(ref.region.1)) 
    message("Reference region 2 set to chr22 by default\n")
  queryC <- unlist(strsplit(x = ref.region.2, split = "[:, -]", 
                            perl = T))
  if (is.na(queryC[2]) | is.na(queryC[3])) {
    subsetC <- GenomicRanges::GRanges(seqnames = queryC[1], 
                                      ranges = IRanges::IRanges(start = 1, end = 3e+08))
  }
  else {
    subsetC <- GenomicRanges::GRanges(seqnames = queryC[1], 
                                      ranges = IRanges::IRanges(start = as.numeric(queryC[2]), 
                                                                end = as.numeric(queryC[3])))
  }
  if (!quiet) 
    message(paste0("Computing coverages in reference regions ", 
                   ref.region.1, ", ", ref.region.2, " and target region ", 
                   target.region, "\n"))
  if (missing(files)) {
    stop("A vector of .bam files paths, a single .bam file path or a folder path containing .bam files must be provided\n")
  }
  else {
    if (length(files) == 1) {
      if (!file.exists(files)) {
        stop("The given path is neither a file or a folder\n")
      }
      else {
        if (file.info(files)$isdir) {
          allfiles <- grep(pattern = ".bam$", x = list.files(files, 
                                                             recursive = T, full.names = T), perl = T, 
                           value = T)
          if (length(allfiles) == 0) 
            stop("There are no files in the given folder\n")
        }
        else {
          if (length(grep(".bam$", files, perl = T, value = T)) == 
              0) 
            stop("The given file is not a .bam file\n")
          allfiles <- files
        }
      }
    }
    else {
      if (any(!file.exists(files))) {
        if (!quiet) 
          message(paste0("The file ", files[!file.exists(files)], 
                         " does not exist\n"))
        files <- files[file.exists(files)]
      }
      if (!(length(grep(".bam$", files, perl = T, value = T)) == 
            length(files))) {
        if (!quiet) 
          message(paste0("The file ", grep(".bam$", files, 
                                           perl = T, value = T, invert = T), " is not a .bam file\n"))
        files <- grep(".bam$", files, perl = T, value = T)
      }
      allfiles <- files
    }
    if (!quiet) 
      message(paste0("Processing ", length(allfiles), " bam file(s)...\n"))
  }
  if (missing(mc.cores)) 
    mc.cores <- 1
  if (mc.cores > length(allfiles)) {
    mc.cores <- length(allfiles)
    if (!quiet) 
      message(paste0("There are more cores than files to be processed. Parameter 'mc.cores' set to ", 
                     mc.cores, "\n"))
  }
  if (missing(exomeTargets)) {
    stop("A path to a exome targets .bed file related to the technology used must be provided\n")
  }
  else {
    if (!file.exists(exomeTargets)) 
      stop("The exome targets .bed file does not exists in the given path.\n")
    if (!(length(grep(".bed$", exomeTargets, perl = T, value = T)) == 
          1)) 
      stop("The exome targets file must be a .bed file\n")
    targets <- read.targets(exomeTargets, skip = skip, quiet = quiet)
    if (GenomeInfoDb::seqlevelsStyle(targets) != "UCSC") {
      GenomeInfoDb::seqlevelsStyle(targets) <- "UCSC"
    }
  }
  targetCov <- parallel::mclapply(X = allfiles, FUN = processBam, 
                                  targets = targets, subset = subsetA, mc.cores = mc.cores)
  ref.1.Cov <- parallel::mclapply(X = allfiles, FUN = processBam, 
                                  targets = targets, subset = subsetB, mc.cores = mc.cores)
  ref.2.Cov <- parallel::mclapply(X = allfiles, FUN = processBam, 
                                  targets = targets, subset = subsetC, mc.cores = mc.cores)
  names(targetCov) <- names(ref.1.Cov) <- names(ref.2.Cov) <- basename(allfiles)
  par <- list(targets = targets, target.region = subsetA, ref.region.1 = subsetB, 
              ref.region.2 = subsetC, files = basename(allfiles))
  coverages <- list(target = targetCov, reference.1 = ref.1.Cov, 
                    reference.2 = ref.2.Cov, par = par)
  class(coverages) <- "MADseqLOY"
  return(coverages)
}