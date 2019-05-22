#' Check the B deviation values in PAR1 and PAR2 regions to check possible Loss of Y events in a folder or files.
#' 
#' checkBdev checks the B deviation values of all the MAD files specified or in
#' a path to check for possible Loss of Y events. The Bdev is computed by default for the 
#' PAR1 and PAR2 regions. The minimum cellularity detected in losses is 26% (bdev.threshold=0.06) by default. Lower values of bdev.threshold can increase
#' the detection of false positives.
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
#' @param top Superior treshold to consider an heterozygous allele. By default is set to 0.85.
#' @param bot Inferior treshold to consider an heterozygous allele. By default is set to 0.15.
#' @param trim trim the fraction (0 to 0.5) of probes to be trimmed when summaryzing LRR.
#' @param mc.cores The number of cores used to perform the function. By default 
#'   is set to 1.
#' @param quiet Should the function not inform about the status of the process. 
#'   By default is FALSE.
#' @param hg Human genome build version.
#' @param pval.sig p-value treshold to be used in the classification test. By default is set to 0.05.
#' @param bdev.threshold bdev threshold to determine if there is BAF split significant enough to call an altered region. By default is set to 0.05.
#' @param ... Other parameters.
#' @return A MADloyBdev object that contains the Bdev values for the two PAR regions for all the files 
#'   analyzed.
#' @export
#' @examples
#' \dontrun{
#' checkBdev(filepath, mc.cores=2)}
checkBdev <- function(object, rsCol = 1, ChrCol = 2, PosCol = 3, LRRCol = 4, BAFCol = 5, 
    top = 0.85, bot = 0.15, trim = 0.1, mc.cores, quiet = FALSE, hg = "hg18", pval.sig = 0.05, 
    bdev.threshold = 0.06, ...) {
    
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
    lrr2ploidy <- function(x) 2 * exp(3 * x/2)
    ploidy2lrr <- function(x) 2 * log(x/2)/3
    
    # Check input-----------------------------------------------------------------
    
    if (missing(object)) {
        stop("A LOY object, a single file path (APT platform and MAD platform), a vector of files paths (MAD platform) or a MAD rawData folder path containing files ready to be processed with MAD (MAD platform) must be provided")
    } else {
        if (inherits(object, "LOY")) {
            if (!quiet) 
                message("Processing the files in the LOY object")
            if (sum(object$par$QCremoved)>=1) {
              if (!quiet)
                message(paste0(sum(object$par$QCremoved), " samples removed due to large LRR standard deviation values"))
              allfiles <- file.path(object$par$path, object$par$files)[!object$par$QCremoved]
              n <- length(object$par$files[!object$par$QCremoved])
              cl <- data.frame(orig = object$res$MADloy[!object$par$QCremoved])
              rownames(cl) <- tools::file_path_sans_ext(object$par$files[!object$par$QCremoved])  
            } else {
              allfiles <- file.path(object$par$path, object$par$files)
              n <- length(object$par$files)
              cl <- data.frame(orig = object$res$MADloy)
              rownames(cl) <- tools::file_path_sans_ext(object$par$files)  
            }
            # process PAR regions -----------------------------------------
            regions <- object$par$regions
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
                cl <- data.frame(orig = rep("Normal", length(allfiles)))
                rownames(cl) <- basename(allfiles)
                # process PAR regions -----------------------------------------
                regions <- fread(system.file("extdata", "references", paste0(hg, 
                  ".par.regions"), package = "MADloy"), header = T, skip = 1, colClasses = c("character", 
                  "character", "numeric", "numeric"), showProgress = FALSE)
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
        ChrCol = ChrCol, PosCol = PosCol, LRRCol = LRRCol, BAFCol = BAFCol, regions = regions, 
        mc.cores = mc.cores, top = top, bot = bot, trim = trim)
    names(data) <- basename(allfiles)
    par <- list(files = basename(allfiles), path = dirname(allfiles), cols = c(rsCol, 
        ChrCol, PosCol, LRRCol, BAFCol), top = top, bot = bot, hg = hg)
    
    if( any(sapply(data, "[[", "check") == "XTR")) {
      fsX <- data.frame(t(sapply(data, "[[", "fsX")), stringsAsFactors = FALSE)
      fsX$Pl <- as.numeric(fsX$Pl)
      XTR <- data.frame(t(sapply(data, "[[", "XTR")), stringsAsFactors = FALSE)
      XTR$Pl <- as.numeric(XTR$Pl)
      
      cl$hetfsX <- round(as.numeric(fsX$n) / as.numeric(fsX$N), 2)
      cl$hetXTR <- round(as.numeric(XTR$n) / as.numeric(XTR$N), 2)
      cl$LRRfsX <- ploidy2lrr(fsX$Pl)
      cl$LRRXTR <- ploidy2lrr(XTR$Pl)
      cl$BdevfsX <- round(unlist(fsX$Bdev), 3)
      cl$BdevXTR <- round(unlist(XTR$Bdev), 3)
      cl$classBdev <- ifelse(XTR$Bdev > bdev.threshold & XTR$Pl < 2, "LOY", "normal")
      cl$classBdev[XTR$Bdev > bdev.threshold & XTR$Pl > 2 ] <- "XYY"
      cl$classBdev[fsX$Bdev > bdev.threshold & fsX$Pl > 2 ] <- "XXY"
      cl$classBdev[cl$hetfsX >= 0.10 ] <- "other"
      cl$class[cl$orig == cl$classBdev] <- cl$classBdev[cl$orig == cl$classBdev]
      cl$class[cl$orig != cl$classBdev] <- "discordant"
      cl$class[cl$classBdev == "other"] <- "other"
      
      cl$LRRCellXTR <- abs(2-unlist(XTR$Pl))*100
      cl$BdevCellXTR <- rep(0, nrow(cl))
      cl$BdevCellXTR[cl$class == "LOY" & !is.na(cl$class)] <- round((2*cl$BdevXTR[cl$class == "LOY" & !is.na(cl$class)])/(0.5+cl$BdevXTR[cl$class == "LOY" & !is.na(cl$class)])*100, 2)
      cl$BdevCellXTR[cl$class == "XYY" & !is.na(cl$class)] <- round((2*cl$BdevXTR[cl$class == "XYY" & !is.na(cl$class)])/(0.5-cl$BdevXTR[cl$class == "XYY" & !is.na(cl$class)])*100, 2)
      Bdev <- list(class = cl$class, data = cl, Bdev = data, par = par)
    } else {
      fsX <- data.frame(t(sapply(data, "[[", "fsX")), stringsAsFactors = FALSE)
      fsX$Pl <- as.numeric(fsX$Pl)
      PAR1 <- data.frame(t(sapply(data, "[[", "PAR1")), stringsAsFactors = FALSE)
      PAR1$Pl <- as.numeric(PAR1$Pl)
      PAR2 <- data.frame(t(sapply(data, "[[", "PAR2")), stringsAsFactors = FALSE)
      PAR2$Pl <- as.numeric(PAR2$Pl)
      XTR <- data.frame(t(sapply(data, "[[", "XTR")), stringsAsFactors = FALSE)
      XTR$Pl <- as.numeric(XTR$Pl)
      
      PARstat <- data.frame(t(sapply(data, function(x) {
        t.test2(x$PAR1$Pl, x$PAR2$Pl, x$PAR1$Plsd, x$PAR2$Plsd, x$PAR1$n, x$PAR2$n, equal.variance = FALSE)
      })))
      cl$LRRfsX <- ploidy2lrr(fsX$Pl)
      cl$LRRPAR1 <- ploidy2lrr(PAR1$Pl)
      cl$LRRPAR2 <- ploidy2lrr(PAR2$Pl)
      cl$LRRXTR <- ploidy2lrr(XTR$Pl)
      
      cl$HetfsX <- round(as.numeric(fsX$n) / as.numeric(fsX$N), 2)
      cl$BdevfsX <- round(unlist(fsX$Bdev), 3)
      cl$BdevPAR1 <- round(unlist(PAR1$Bdev), 3)
      cl$BdevPAR2 <- round(unlist(PAR2$Bdev), 3)
      cl$BdevXTR <- round(unlist(XTR$Bdev), 3)
      
      cl$classBdev <- ifelse(PAR1$Bdev > bdev.threshold & PAR2$Bdev > bdev.threshold & XTR$Bdev > bdev.threshold & PAR1$Pl < 2 & PAR2$Pl < 2 , "LOY", "normal")
      cl$classBdev[PAR1$Bdev > bdev.threshold & PAR2$Bdev > bdev.threshold & XTR$Bdev > bdev.threshold & PAR1$Pl > 2 & PAR2$Pl > 2 ] <- "XYY"
      cl$classBdev[cl$HetfsX >= 0.10] <- "other"
      cl$class[cl$orig == cl$classBdev] <- cl$classBdev[cl$orig == cl$classBdev]
      cl$class[cl$orig != cl$classBdev] <- "discordant"
      cl$class[cl$classBdev == "other"] <- "other"
      cl$LRRCellPAR1 <- abs(2-unlist(PAR1$Pl))*100
      cl$LRRCellPAR2 <- abs(2-unlist(PAR2$Pl))*100
      cl$LRRCellXTR <- abs(2-unlist(XTR$Pl))*100
      cl$BdevCellXTR <- cl$BdevCellPAR2 <- cl$BdevCellPAR1 <- rep(0, nrow(cl))
      cl$BdevCellPAR1[cl$class == "LOY" & !is.na(cl$class)] <- round((2*cl$BdevPAR1[cl$class == "LOY" & !is.na(cl$class)])/(0.5+cl$BdevPAR1[cl$class == "LOY" & !is.na(cl$class)])*100, 2)
      cl$BdevCellPAR2[cl$class == "LOY" & !is.na(cl$class)] <- round((2*cl$BdevPAR2[cl$class == "LOY" & !is.na(cl$class)])/(0.5+cl$BdevPAR2[cl$class == "LOY" & !is.na(cl$class)])*100, 2)
      cl$BdevCellXTR[cl$class == "LOY" & !is.na(cl$class)] <- round((2*cl$BdevXTR[cl$class == "LOY" & !is.na(cl$class)])/(0.5+cl$BdevXTR[cl$class == "LOY" & !is.na(cl$class)])*100, 2)
      cl$BdevCellPAR1[cl$class == "XYY" & !is.na(cl$class)] <- round((2*cl$BdevPAR1[cl$class == "XYY" & !is.na(cl$class)])/(0.5-cl$BdevPAR1[cl$class == "XYY" & !is.na(cl$class)])*100, 2)
      cl$BdevCellPAR2[cl$class == "XYY" & !is.na(cl$class)] <- round((2*cl$BdevPAR2[cl$class == "XYY" & !is.na(cl$class)])/(0.5-cl$BdevPAR2[cl$class == "XYY" & !is.na(cl$class)])*100, 2)
      cl$BdevCellXTR[cl$class == "XYY" & !is.na(cl$class)] <- round((2*cl$BdevXTR[cl$class == "XYY" & !is.na(cl$class)])/(0.5-cl$BdevXTR[cl$class == "XYY" & !is.na(cl$class)])*100, 2)
      cl$balanced <- ifelse(PARstat$p.value > pval.sig * 10/nrow(PARstat), "balancedPAR", 
                            "unbalancedPAR")
      cl$balanced[cl$orig == "LOY" & PARstat$p.value < pval.sig/n & PAR1$Pl > PAR2$Pl] <- "LOYq"
      cl$balanced[cl$orig == "LOY" & PARstat$p.value < pval.sig/n & PAR1$Pl < PAR2$Pl] <- "LOYp"
      cl$balanced[cl$orig == "XYY" & PARstat$p.value < pval.sig/n & PAR1$Pl < PAR2$Pl] <- "XYYp"
      cl$balanced[cl$orig == "XYY" & PARstat$p.value < pval.sig/n & PAR1$Pl < PAR2$Pl] <- "XYYq"
      Bdev <- list(class = cl$class, data = cl, prob = PARstat, Bdev = data, par = par)
    }

# Merge both information (JRG - 22nd May)        
    dat1 <- data.frame(MADthres = object$res$MADloy,
                       Fosberg = object$res$Fosberg, 
                       mLRR = object$res$continous)
    rownames(dat1) <- gsub(".txt", "", rownames(object$res))
    dat2 <- data.frame(orig=Bdev$data$orig, class=Bdev$data$class)
    rownames(dat2) <- gsub(".txt", "", rownames(Bdev$data))
    
    ans <- merge(dat1, dat2, by="row.names", all.x=TRUE)
    ans$MADloy <- ans$class
    ans$MADloy[is.na(ans$class)] <- ans$MADthres[is.na(ans$class)]
    ans$MADloy <- relevel(ans$MADloy, ref="normal")
    rownames(ans) <- ans$Row.names
    names(ans)[1] <- "id"
    
    out <- list(calling=ans, Bdev=Bdev)
    class(out) <- "MADloyBdev"
    return(out)
}
