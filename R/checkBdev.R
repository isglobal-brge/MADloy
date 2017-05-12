#' Check the B deviation values in PAR1 and PAR2 regions to check possible Loss of Y events in a folder or files.
#' 
#' checkBdev checks the B deviation values of all the MAD files specified or in
#' a path to check for possible Loss of Y events. The Bdev is computed by default for the 
#' PAR1 and PAR2 regions.
#' 
#' @seealso \code{\link{getLOY}} to process results from \code{MADloy}
#' @param files A LOY class object to check the LOY and uncertain samples or a single file path (APT platform and MAD platform), a vector of 
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
checkBdev <- function( object, rsCol = 1, ChrCol = 2, PosCol = 3, LRRCol= 4, BAFCol = 5, top = 0.9, bot = 0.1, trim = 0, mc.cores, quiet = FALSE, hg="hg18", pval.sig = 0.05, ...) {
  
  #two-sample t-test from https://stats.stackexchange.com/questions/30394/how-to-perform-two-sample-t-tests-in-r-by-inputting-sample-statistics-rather-tha
  t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
  {
    if( equal.variance==FALSE ) 
    {
      se <- sqrt( (s1^2/n1) + (s2^2/n2) )
      # welch-satterthwaite df
      df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    } else
    {
      # pooled standard deviation, scaled by the sample sizes
      se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
      df <- n1+n2-2
    }      
    t <- (m1-m2-m0)/se 
    dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
    names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
    return(dat) 
  }
  
  # process PAR regions -----------------------------------------
  
  par <- fread(system.file("data", paste0(hg, ".par.regions"), package = "MADloy"), header=T, skip=1, colClasses = c("character", "character", "numeric", "numeric"), showProgress = FALSE)
  
  subset <- GenomicRanges::GRanges(seqnames = gsub("chr", "", par[par$chromosome == "Y"]$chromosome), ranges = IRanges::IRanges(start = par[par$chromosome == "Y"]$start, end = par[par$chromosome == "Y"]$end))
  
  # Check input-----------------------------------------------------------------
  
  if (missing(object)) {
    stop("A LOY object, a single file path (APT platform and MAD platform), a vector of files paths (MAD platform) or a MAD rawData folder path containing files ready to be processed with MAD (MAD platform) must be provided")
  } else {
    if (inherits(object, "LOY")) {
      if (!quiet) message("Processing the files in the LOY object")
      allfiles <- file.path(object$par$path, object$par$files)[object$prob <= 0.05/length(object$par$files)]
    } else {
      if (length(object) == 1) {
        if (!file.exists(object)) {
          stop("The given object is neither a file path or a folder")
        } else {
          if (file.info(object)$isdir) {
            allfiles <- list.files(object, recursive = T, full.names = T)
            if (length(allfiles) == 0) 
              stop("There are no files in the given folder")
          } else {
            allfiles <- object
          }
        }
      } else {
        if (any(!file.exists(object))) {
          if (!quiet) 
            message(paste0("The file ", object[!file.exists(object)], " does not exist"))
          object <- object[file.exists(object)]
        }
        allfiles <- object
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
  
  # Get Bdev summary ------------------------------------------------------------------
  
  data <- parallel::mclapply(X = allfiles, FUN = MADloy:::processBdevMAD, rsCol = rsCol, ChrCol = ChrCol, 
    PosCol = PosCol, LRRCol = LRRCol, BAFCol = BAFCol, query = subset, mc.cores = mc.cores, top = top, bot = bot, trim = 0.1 )
  names(data) <- basename(allfiles)
  par <- list(files = basename(allfiles),
              path = dirname(allfiles), cols = c(rsCol, ChrCol, PosCol, LRRCol, BAFCol),
              top = top,
              bot = bot,
              hg = hg)
  
  p <- data.frame(t(sapply(data, "[[", "p")), stringsAsFactors = FALSE)
  p$Pl <- as.numeric(p$Pl)
  q <- data.frame(t(sapply(data, "[[", "q")), stringsAsFactors = FALSE)
  q$Pl <- as.numeric(q$Pl)
  pqstat <- data.frame(t(sapply(data, function(x){t.test2(x$p$Pl, x$q$Pl, x$p$Plsd, x$q$Plsd, x$p$n, x$q$n, equal.variance = FALSE)})))
  
  cl <- data.frame(orig = object$class[object$prob <= 0.05/length(object$par$files)])
  
  cl$pq <- ifelse(pqstat$p.value > pval.sig*10/nrow(pqstat), "balancedpq", "unbalancedpq")
  cl$pqn <- ifelse(pqstat$p.value > pval.sig, "balancedpq", "unbalancedpq")
  cl$pq[cl$orig=="LOY" &  pqstat$p.value < pval.sig/length(object$par$files) & p$Pl > q$Pl] <- "LOYq"
  cl$pq[cl$orig=="LOY" &  pqstat$p.value < pval.sig/length(object$par$files) & p$Pl < q$Pl] <- "LOYp"
  cl$pq[cl$orig=="XYY" &  pqstat$p.value < pval.sig/length(object$par$files) & p$Pl < q$Pl] <- "XYYp"
  cl$pq[cl$orig=="XYY" &  pqstat$p.value < pval.sig/length(object$par$files) & p$Pl < q$Pl] <- "XYYq"
  

  
  Bdev <- list(cl = cl, Bdev = data, par = par)
  
  class(Bdev) <- "MADloyBdev"
  return(Bdev)
} 
