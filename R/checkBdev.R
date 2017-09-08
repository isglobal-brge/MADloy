#' Check the B deviation values in PAR1 and PAR2 regions to check possible Loss of Y events in a folder or files.
#' 
#' checkBdev checks the B deviation values of all the MAD files specified or in
#' a path to check for possible Loss of Y events. The Bdev is computed by default for the 
#' PAR1 and PAR2 regions.
#' 
#' @seealso \code{\link{getLOY}} to process results from \code{MADloy}
#' @param object A LOY class object to check the LOY and uncertain samples or a single file path (APT platform and MAD platform), a vector of 
#'   file paths (MAD platform) or a MAD rawData folder path containing files 
#'   ready to be processed with MAD (MAD platform).
#' @param rsCol The position of the column with the SNP identifier.
#' @param ChrCol The position of the column with the Chromosome field.
#' @param PosCol The position of the column with the Position field.
#' @param LRRCol The position of the column with the LRR field.
#' @param BAFCol The position of the column with the BAF field.
#' @param top Superior treshold to consider an heterozygous allele.
#' @param bot Inferior treshold to consider an heterozygous allele.
#' @param trim trim the fraction (0 to 0.5) of probes to be trimmed when summaryzing LRR.
#' @param mc.cores The number of cores used to perform the function. By default 
#'   is set to 1.
#' @param quiet Should the function not inform about the status of the process. 
#'   By default is FALSE.
#' @param hg Human genome build version.
#' @param pval.sig p-value treshold to be used in the classification test.
#' @param ... Other parameters.
#' @return A MADloyBdev object that contains the Bdev values for the two PAR regions for all the files 
#'   analyzed.
#' @export
#' @examples
#' \dontrun{
#' checkBdev(filepath, mc.cores=2)}
checkBdev <- function(object, rsCol = 1, ChrCol = 2, PosCol = 3, LRRCol = 4, BAFCol = 5, 
    top = 0.7, bot = 0.3, trim = 0.1, mc.cores, quiet = FALSE, hg = "hg18", pval.sig = 0.05, 
    ...) {
    
    # two-sample t-test from
    # https://stats.stackexchange.com/questions/30394/how-to-perform-two-sample-t-tests-in-r-by-inputting-sample-statistics-rather-tha
    t.test2 <- function(m1, m2, s1, s2, n1, n2, m0 = 0, equal.variance = FALSE) {
        if (equal.variance == FALSE) {
            se <- sqrt((s1^2/n1) + (s2^2/n2))
            # welch-satterthwaite df
            df <- ((s1^2/n1 + s2^2/n2)^2)/((s1^2/n1)^2/(n1 - 1) + (s2^2/n2)^2/(n2 - 
                1))
        } else {
            # pooled standard deviation, scaled by the sample sizes
            se <- sqrt((1/n1 + 1/n2) * ((n1 - 1) * s1^2 + (n2 - 1) * s2^2)/(n1 + 
                n2 - 2))
            df <- n1 + n2 - 2
        }
        t <- (m1 - m2 - m0)/se
        dat <- c(m1 - m2, se, t, 2 * stats::pt(-abs(t), df))
        names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
        return(dat)
    }
    
    
    # Check input-----------------------------------------------------------------
    
    if (missing(object)) {
        stop("A LOY object, a single file path (APT platform and MAD platform), a vector of files paths (MAD platform) or a MAD rawData folder path containing files ready to be processed with MAD (MAD platform) must be provided")
    } else {
        if (inherits(object, "LOY")) {
            if (!quiet) 
                message("Processing the files in the LOY object")
            allfiles <- file.path(object$par$path, object$par$files)[object$prob <= 
                0.05/length(object$par$files)]
            n <- length(object$par$files)
            cl <- data.frame(orig = object$class[object$prob <= 0.05/n])
            # process PAR regions -----------------------------------------
            subset <- GenomicRanges::GRanges(seqnames = gsub("chr", "", object$par$regions[object$par$regions$chromosome == 
                "Y"]$chromosome), ranges = IRanges::IRanges(start = object$par$regions[object$par$regions$chromosome == 
                "Y"]$start, end = object$par$regions[object$par$regions$chromosome == 
                "Y"]$end))
        } else {
            if (length(object) == 1) {
                if (!file.exists(object)) {
                  stop("The given object is neither a file path or a folder")
                } else {
                  if (file.info(object)$isdir) {
                    allfiles <- list.files(object, recursive = T, full.names = T)
                    n <- length(allfiles)
                    cl <- data.frame(orig = allfiles)
                    if (length(allfiles) == 0) 
                      stop("There are no files in the given folder")
                  } else {
                    allfiles <- object
                    n <- length(allfiles)
                    cl <- data.frame(orig = allfiles)
                  }
                }
            } else {
                if (any(!file.exists(object))) {
                  if (!quiet) 
                    message(paste0("The file ", object[!file.exists(object)], " does not exist"))
                  object <- object[file.exists(object)]
                }
                allfiles <- object
                n <- length(allfiles)
                cl <- data.frame(getLOY = allfiles)
                # process PAR regions -----------------------------------------
                regions <- fread(system.file("extdata", "references", paste0(hg, 
                  ".par.regions"), package = "MADloy"), header = T, skip = 1, colClasses = c("character", 
                  "character", "numeric", "numeric"), showProgress = FALSE)
                subset <- GenomicRanges::GRanges(seqnames = gsub("chr", "", regions[regions$chromosome == "Y"]$chromosome), 
                                                 ranges = IRanges::IRanges(start = regions[regions$chromosome == "Y"]$start,
                                                                           end = regions[regions$chromosome == "Y"]$end))
            }
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
    
    # Get Bdev summary
    # ------------------------------------------------------------------
    
    data <- parallel::mclapply(X = allfiles, FUN = processBdevMAD, rsCol = rsCol, 
        ChrCol = ChrCol, PosCol = PosCol, LRRCol = LRRCol, BAFCol = BAFCol, query = subset, 
        mc.cores = mc.cores, top = top, bot = bot, trim = trim)
    names(data) <- basename(allfiles)
    par <- list(files = basename(allfiles), path = dirname(allfiles), cols = c(rsCol, 
        ChrCol, PosCol, LRRCol, BAFCol), top = top, bot = bot, hg = hg)
    
    p <- data.frame(t(sapply(data, "[[", "p")), stringsAsFactors = FALSE)
    p$Pl <- as.numeric(p$Pl)
    q <- data.frame(t(sapply(data, "[[", "q")), stringsAsFactors = FALSE)
    q$Pl <- as.numeric(q$Pl)
    pqstat <- data.frame(t(sapply(data, function(x) {
        t.test2(x$p$Pl, x$q$Pl, x$p$Plsd, x$q$Plsd, x$p$n, x$q$n, equal.variance = FALSE)
    })))
    cl$adjusted_p <- ifelse(pqstat$p.value > pval.sig * 10/nrow(pqstat), "balancedpq", 
        "unbalancedpq")
    cl$adjusted_p[cl$getLOY == "LOY" & pqstat$p.value < pval.sig/n & p$Pl > q$Pl] <- "LOYq"
    cl$adjusted_p[cl$getLOY == "LOY" & pqstat$p.value < pval.sig/n & p$Pl < q$Pl] <- "LOYp"
    cl$adjusted_p[cl$getLOY == "XYY" & pqstat$p.value < pval.sig/n & p$Pl < q$Pl] <- "XYYp"
    cl$adjusted_p[cl$getLOY == "XYY" & pqstat$p.value < pval.sig/n & p$Pl < q$Pl] <- "XYYq"
    Bdev <- list(class = cl, prob = pqstat, Bdev = data, par = par)
    
    class(Bdev) <- "MADloyBdev"
    return(Bdev)
}
