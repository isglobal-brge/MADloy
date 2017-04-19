#' Check the B deviation values in PAR1 and PAR2 regions to check possible Loss of Y events in a folder or files.
#' 
#' checkBdev checks the B deviation values of all the MAD files specified or in
#' a path to check for possible Loss of Y events. The Bdev is computed by default for the 
#' PAR1 and PAR2 regions.
#' 
#' @seealso \code{\link{getLOY}} to process results from \code{MADloy}
#' @param files A single file path (APT platform and MAD platform), a vector of 
#'   file paths (MAD platform) or a MAD rawData folder path containing files 
#'   ready to be processed with MAD (MAD platform).
#' @param rsCol The position of the column with the SNP identifier.
#' @param ChrCol The position of the column with the Chromosome field.
#' @param PosCol The position of the column with the Position field.
#' @param BAFCol The position of the column with the BAF identifier.
#' @param top Superior treshold to consider an heterozygous allele.
#' @param bot Inferior treshold to consider an heterozygous allele.
#' @param mc.cores The number of cores used to perform the function. By default 
#'   is set to 1.
#' @param quiet Should the function not inform about the status of the process. 
#'   By default is FALSE.
#' @param ... Other parameters
#' @return A MADloyBdev object that contains the Bdev values for the two PAR regions for all the files 
#'   analyzed.
#' @export
#' @examples
#' \dontrun{
#' checkBdev(filepath, mc.cores=2)}
checkBdev <- function(files, rsCol = 1, ChrCol = 2, PosCol = 3, BAFCol = 4, top=0.9, bot=0.1, mc.cores, quiet = FALSE, hg="hg18", ...) {
  
  # process PAR regions -----------------------------------------
  
  par <- fread(system.file("data", paste0(hg, ".par.regions"), package = "MADloy"), header=T, skip=1, colClasses = c("character", "character", "numeric", "numeric"), showProgress = FALSE)
  
  subset <- GenomicRanges::GRanges(seqnames = gsub("chr", "", par[par$chromosome == "Y"]$chromosome), ranges = IRanges::IRanges(start = par[par$chromosome == "Y"]$start, end = par[par$chromosome == "Y"]$end))
  
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
  
  # Get Bdev summary ------------------------------------------------------------------
  
  Bdev <- parallel::mclapply(X = allfiles, FUN = processBdevMAD, rsCol = rsCol, ChrCol = ChrCol, 
    PosCol = PosCol, BAFCol = BAFCol, query = subset, mc.cores = mc.cores, top = top, bot = bot)
  names(Bdev) <- basename(allfiles)
  par <- list(files = basename(allfiles),
              path = dirname(allfiles), cols = c(rsCol, ChrCol, PosCol, LRRCol, BAFCol),
              top = top,
              bot = bot,
              hg = hg)
  Bdev <- list(Bdev = Bdev, par = par)
  class(Bdev) <- "MADloyBdev"
  return(Bdev)
} 
