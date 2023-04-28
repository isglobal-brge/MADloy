#' Check sex for SNP array samples
#' 
#' @param files A single file path (APT platform and MAD platform), a vector of 
#'   file paths (MAD platform) or a MAD rawData folder path containing files 
#'   ready to be processed with MAD (MAD platform).
#' @param rsCol The position of the column with the SNP identifier.
#' @param ChrCol The position of the column with the Chromosome field.
#' @param PosCol The position of the column with the Position field.
#' @param LRRCol The position of the column with the LRR identifier.
#' @param mc.cores The number of cores used to perform the function. By default 
#'   is set to 1.
#' @param trim the fraction (0 to 0.5) of probes to be trimmed when summaryzing LRR. 
#' This argument tries to control the effect of having CNVs across genome. Default is 0.2.
#' @param quiet Should the function not inform about the status of the process. 
#'   By default is FALSE.
#' @return A checkSex object that contains the trimmed mean LRR in chrX and chrY in
#'   the analyzed samples, and its classification as males and females in the class variableinside the object.
#' @export
#' @examples
#' \dontrun{
#' checkSex(filepath, mc.cores=20)}
checkSex <- function(files, rsCol = 1, ChrCol = 2, PosCol = 3, LRRCol = 4, mc.cores, 
    trim = 0.2, quiet = FALSE) {
    
    lrr2ploidy <- function(x) 2 * exp(3 * x/2)
    ploidy2lrr <- function(x) 2 * log(x/2)/3
    ploidy <- c(0.001, 1, 2, 3, 4, 5)
    lrrploidy <- ploidy2lrr(ploidy)
    getXY <- function(x, rsCol, ChrCol, PosCol, LRRCol, trim) {
        dat <- data.table::fread(x, showProgress = FALSE, sep = "\t")
        data.table::setnames(dat, colnames(dat[, c(rsCol, ChrCol, PosCol, LRRCol), 
            with = F]), c("Name", "Chr", "Position", "Log.R.Ratio"))
        XYsummary <- list()
        XYsummary$X <- mean(dat$Log.R.Ratio[dat$Chr == "X"], na.rm = T, trim = trim)
        XYsummary$Y <- mean(dat$Log.R.Ratio[dat$Chr == "Y"], na.rm = T, trim = trim)
        return(XYsummary)
    }
    
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
    
    res <- parallel::mclapply(X = allfiles, FUN = getXY, rsCol = rsCol, ChrCol = ChrCol,
        PosCol = PosCol, LRRCol = LRRCol, mc.cores = mc.cores, trim = trim)
    
    data <- do.call(rbind, res)
    rownames(data) <- basename(allfiles)
    data <- as.data.frame(apply(data, 2, unlist))
    omitted <- data[which(!complete.cases(data)),]
    data <- data[which(complete.cases(data)),]

    ## Clustering
    
    centers <- rbind(c(lrrploidy[3], lrrploidy[1]), c(lrrploidy[2], lrrploidy[2]))
    kmeansres <- try(stats::kmeans(lrr2ploidy(data), 
                                   centers = lrr2ploidy(centers)), TRUE)
    if (!inherits(kmeansres, "try-error")) {
      offsetY <- stats::median(data[kmeansres$cluster == 2, ]$Y) - lrrploidy[2]
      offsetX <- stats::median(data[kmeansres$cluster == 1, ]$X) - lrrploidy[3]
    
      tmp <- data
      tmp$X <- ploidy2lrr(lrr2ploidy(tmp$X) + (2 - lrr2ploidy(offsetX)))
      tmp$Y <- tmp$Y - offsetY
      
      if(sum(!complete.cases(tmp))>0){
          omitted <- rbind(omitted, tmp[which(!complete.cases(tmp)),])
          data <- data[which(complete.cases(tmp)),]
          tmp <- tmp[which(complete.cases(tmp)),]
      }
    
      kmeansres <- stats::kmeans(lrr2ploidy(tmp), centers = lrr2ploidy(centers))
      class <- as.factor(ifelse(kmeansres$cluster == 1, "FEMALE", "MALE"))
    }
    else {
      offsetY <- stats::median(data$Y) - lrrploidy[2]
      offsetX <- 0
      class <- factor(rep("MALE", length(data$Y)), 
                      levels=c("FEMALE", "MALE"))
      warning("All samples have the same gender status.")
    }
    
    par <- list()
    par$trim <- trim
    par$offsetY <- offsetY
    par$offsetX <- offsetX
    par$path <- files
    par$omitted <- rownames(omitted)
    par$files <- allfiles
    
    res <- list(data = data,
                class = class, par = par)
    class(res) <- "checkSex"
    return(res)
}
