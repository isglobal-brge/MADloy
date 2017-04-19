#' Check the mean LRR values in msY region to detect Loss of Y events in a folder or files.
#' 
#' madloy check the median log R ratio(LRR) of all the MAD files specified or in
#' a path to detect Loss of Y events. The LRR is computed by default for the 
#' autosomes and msY chrY:2694521-59034049 (hg19/GRCh37).
#' 
#' @seealso \code{\link{getLOY}} to process results from \code{MADloy}
#' @param files A single file path (APT platform and MAD platform), a vector of 
#'   file paths (MAD platform) or a MAD rawData folder path containing files 
#'   ready to be processed with MAD (MAD platform).
#' @param target.region The chromosome or region to be compared with the other 
#'   regions. By default is the region chrY:2694521-59034049 (hg19/GRCh37) but it can be 
#'   changed.
#' @param ref.region If declared, the chromosome or region to be compared with the Y 
#'   region in UCSC style (i.e. "chr21" or "chr21:1000-10000").
#' @param rsCol The position of the column with the SNP identifier.
#' @param ChrCol The position of the column with the Chromosome field.
#' @param PosCol The position of the column with the Position field.
#' @param LRRCol The position of the column with the LRR identifier.
#' @param trim trim the fraction (0 to 0.5) of probes to be trimmed when summaryzing LRR.
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
madloy <- function(files, target.region = "chrY:2694521-59034049", ref.region="Autosomes",
  rsCol = 1, ChrCol = 2, PosCol = 3, LRRCol = 4, trim=0, mc.cores, quiet = FALSE, hg="hg19", ...) {
  
  chrSizes <- fread(system.file("data", paste0(hg, ".chrom.sizes"), package = "MADloy"), header=T, skip="#", colClasses = c("character", "numeric"), showProgress = FALSE)
  
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
  
  if (missing(ref.region)) {
    message("Using all autosomes as Reference region")
    subsetB <- GenomicRanges::GRanges(seqnames = chrSizes$chromosome[1:22], ranges = IRanges::IRanges(0, chrSizes$size[1:22]))
  } else {
    queryB <- unlist(strsplit(x = ref.region, split = "[:, -]", perl = T))
    if (is.na(queryB[2]) | is.na(queryB[3])) {
      subsetB <- GenomicRanges::GRanges(seqnames = gsub("chr", "", queryB[1]), ranges = IRanges::IRanges(start = 1, 
                                                                                                        end = 3e+08))
      if (queryB[1] == "Autosomes") subsetB <- GenomicRanges::GRanges(seqnames = chrSizes$chromosome[1:22], ranges = IRanges::IRanges(0, chrSizes$size[1:22]))
    } else {
      subsetB <- GenomicRanges::GRanges(seqnames = gsub("chr", "", queryB[1]), ranges = IRanges::IRanges(start = as.numeric(queryB[2]), 
                                                                                                        end = as.numeric(queryB[3])))
    }
  }
  
  if (!quiet) 
    if (missing(ref.region)) message(paste0("LRR median will be computed in all autosomal chromosomes as reference regions and target region ", target.region, "\n")) else message(paste0("LRR median will be computed in reference region ", ref.region, " and target region ", target.region, "\n"))
  
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
  
  # Get LRR summary ------------------------------------------------------------------
  
  targetLRR <- parallel::mclapply(X = allfiles, FUN = processMAD, rsCol = rsCol, ChrCol = ChrCol, 
    PosCol = PosCol, LRRCol = LRRCol, query = subsetA, mc.cores = mc.cores, trim=trim)
  refLRR <- parallel::mclapply(X = allfiles, FUN = processMAD, rsCol = rsCol, ChrCol = ChrCol, 
    PosCol = PosCol, LRRCol = LRRCol, query = subsetB, mc.cores = mc.cores, trim=trim)
  names(targetLRR) <- names(refLRR) <- basename(allfiles)
  par <- list(target.region = subsetA, ref.region = subsetB, 
              files = basename(allfiles),
              path = dirname(allfiles), cols = c(rsCol, ChrCol, PosCol, LRRCol),
              trim = trim)
  LRRsummary <- list(target = targetLRR, reference = refLRR, par = par)
  class(LRRsummary) <- "MADloy"
  return(LRRsummary)
} 
