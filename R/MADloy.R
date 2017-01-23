#' Check the mean LRR values to detect Loss of Y events in a folder or files.
#' 
#' MADloy check the median log R ratio(LRR) of all the MAD files specified or in
#' a path to detect Loss of Y events. The LRR is computed by default for the 
#' regions chr21, chr22 and chrY:2694521-59034049 (hg19/GRCh37).
#' 
#' @seealso \code{\link{getLOY}} to process results from \code{MADloy}
#' @param files A single file path (APT platform and MAD platform), a vector of 
#'   file paths (MAD platform) or a MAD rawData folder path containing files 
#'   ready to be processed with MAD (MAD platform).
#' @param target.region The chromosome or region to be compared with the other 
#'   regions. By default is the region chrY:2694521-59034049 (hg19/GRCh37) but it can be 
#'   changed.
#' @param ref.region.1 First chromosome or region to be compared with the Y 
#'   region in UCSC style (i.e. "chr21" or "chr21:1000-10000").
#' @param ref.region.2 Second chromosome or region to be compared with the Y 
#'   region in UCSC style (i.e. "chr22" or "chr22:1000-10000").
#' @param rsCol The position of the column with the SNP identifier.
#' @param ChrCol The position of the column with the Chromosome field.
#' @param PosCol The position of the column with the Position field.
#' @param LRRCol The position of the column with the LRR identifier.
#' @param mc.cores The number of cores used to perform the function. By default 
#'   is set to 1.
#' @param quiet Should the function not inform about the status of the process. 
#'   By default is FALSE.
#' @param ... Other parameters
#' @return A MADloy object that contains the LRR means for all the files 
#'   analyzed.
#' @export
#' @examples
#' \dontrun{
#' madloy(filepath, mc.cores=2)}
madloy <- function(files, target.region = "chrY:2694521-59034049", ref.region.1 = "chr21", ref.region.2 = "chr22",
  rsCol = 1, ChrCol = 2, PosCol = 3, LRRCol = 4, mc.cores, quiet = FALSE, ...) {
  
  # Check target and reference regions -----------------------------------------
  if (missing(target.region)) 
    message("Targeted region set to chrY:2694521-59034049 by default")
  queryA <- unlist(strsplit(x = target.region, split = "[:, -]", perl = T))
  if (is.na(queryA[2]) | is.na(queryA[3])) {
    subsetA <- GenomicRanges::GRanges(seqnames = gsub("chr", "", queryA[1]), ranges = IRanges::IRanges(start = 1, 
      end = 3e+08))
  } else {
    subsetA <- GenomicRanges::GRanges(seqnames = gsub("chr", "", queryA[1]), ranges = IRanges::IRanges(start = as.numeric(queryA[2]), 
      end = as.numeric(queryA[3])))
  }
  
  if (missing(ref.region.1)) 
    message("Reference region 1 set to chr21 by default")
  queryB <- unlist(strsplit(x = ref.region.1, split = "[:, -]", perl = T))
  if (is.na(queryB[2]) | is.na(queryB[3])) {
    subsetB <- GenomicRanges::GRanges(seqnames =gsub("chr", "", queryB[1]), ranges = IRanges::IRanges(start = 1, 
      end = 3e+08))
  } else {
    subsetB <- GenomicRanges::GRanges(seqnames =gsub("chr", "", queryB[1]), ranges = IRanges::IRanges(start = as.numeric(queryB[2]), 
      end = as.numeric(queryB[3])))
  }
  
  if (missing(ref.region.2)) 
    message("Reference region 2 set to chr22 by default")
  queryC <- unlist(strsplit(x = ref.region.2, split = "[:, -]", perl = T))
  if (is.na(queryC[2]) | is.na(queryC[3])) {
    subsetC <- GenomicRanges::GRanges(seqnames =gsub("chr", "", queryC[1]), ranges = IRanges::IRanges(start = 1, 
      end = 3e+08))
  } else {
    subsetC <- GenomicRanges::GRanges(seqnames =gsub("chr", "", queryC[1]), ranges = IRanges::IRanges(start = as.numeric(queryC[2]), 
      end = as.numeric(queryC[3])))
  }
  if (!quiet) 
    message(paste0("LRR median will be computed in reference regions ", ref.region.1, 
      ", ", ref.region.2, " and target region ", target.region, "\n"))
  
  
  # Check input-----------------------------------------------------------------
  
  if (missing(files)) {
    stop("A single file path (APT platform and MAD platform), a vector of files paths (MAD platform) or a MAD rawData folder path containing files ready to be processed with MAD (MAD platform) must be provided")
  } else {
    if (length(files) == 1) {
      if (!file.exists(files)) {
        stop("The given path is neither a file or a folder")
      } else {
        if (file.info(files)$isdir) {
          allfiles <- list.files(files, recursive = T, full.names = T)
          if (length(allfiles) == 0) 
          stop("There are no files in the given folder")
        } else {
          allfiles <- files
        }
      }
    } else {
      if (any(!file.exists(files))) {
        if (!quiet) 
          message(paste0("The file ", files[!file.exists(files)], " does not exist"))
        files <- files[file.exists(files)]
      }
      allfiles <- files
    }
    if (!quiet) 
      message(paste0("Processing ", length(allfiles), " file(s)..."))
  }
  
  # Check the number of cores---------------------------------------------------
  
  if (missing(mc.cores)) 
    mc.cores <- 1
  if (mc.cores > length(allfiles)) {
    mc.cores <- length(allfiles)
    if (!quiet) 
      message(paste0("There are more cores than files to be processed. Parameter 'mc.cores' set to ", 
        mc.cores))
  }
  
  # Check LRRs------------------------------------------------------------------
  
  targetLRR <- parallel::mclapply(X = allfiles, FUN = processMAD, rsCol = rsCol, ChrCol = ChrCol, 
    PosCol = PosCol, LRRCol = LRRCol, query = subsetA, mc.cores = mc.cores)
  ref1LRR <- parallel::mclapply(X = allfiles, FUN = processMAD, rsCol = rsCol, ChrCol = ChrCol, 
    PosCol = PosCol, LRRCol = LRRCol, query = subsetB, mc.cores = mc.cores)
  ref2LRR <- parallel::mclapply(X = allfiles, FUN = processMAD, rsCol = rsCol, ChrCol = ChrCol, 
    PosCol = PosCol, LRRCol = LRRCol, query = subsetC, mc.cores = mc.cores)
  names(targetLRR) <- names(ref1LRR) <- names(ref2LRR) <- basename(allfiles)
  par <- list(target.region = subsetA, ref.region.1 = subsetB, 
              ref.region.2 = subsetC, files = basename(allfiles),
              path = dirname(allfiles), cols = c(rsCol, ChrCol, PosCol, LRRCol))
  LRRmedians <- list(target = targetLRR, reference.1 = ref1LRR, reference.2 = ref2LRR, 
    par = par)
  class(LRRmedians) <- "MADloy"
  return(LRRmedians)
} 
