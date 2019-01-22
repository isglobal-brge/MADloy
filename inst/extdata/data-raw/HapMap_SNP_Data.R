# -------------------------------------- HapMap Exome Data

library(data.table)

FTPpath <- "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/"

download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20130502.phase3.exome.alignment.index", destfile="/tmp/20130502.phase3.exome.alignment.index")
HapMapLOYExomeSamples <- fread("/tmp/20130502.phase3.exome.alignment.index")

HapMapLOYExomeSamples <- HapMapLOYExomeSamples[grep(".mapped.ILLUMINA.bwa.CEU.", HapMapLOYExomeSamples$`BAM FILE`, fixed=T)]
exome

for (i in 1:nrow(HapMapLOYsamples)){
  system(paste0("wget ", FTPpath, HapMapLOYsamples$"BAM FILE"[i]))
  system(paste0("wget ", FTPpath, HapMapLOYsamples$"BAI FILE"[i]))
}

gender <- c("Female", "Male", "Female", "Male", "Female", "Male", "Female", "Male", "Male", "Female", "Male", "Female", "Female", "Male", "Male", "Male", "Female", "Male", "Female", "Male", "Female", "Male", "Female", "Male", "Male", "Male", "Female", "Female", "Female", "Male", "Female", "Male", "Female", "Male", "Female", "Male", "Female", "Male", "Female", "Male", "Female", "Male", "Female", "Male", "Female")
exomeMaleSamples <- sapply(strsplit(HapMapLOYsamples[gender == "Male"]$"BAM FILE", "/", fixed=T), "[[", 2)

files <- list.files("/NewLacie_CRW10082/DATASETS/STUDY/Hapmap/SNP/rawData", full.names=T)
SNPsamples <- basename(files)
SNPsamples <- gsub(".txt", "", SNPsamples)

files[SNPsamples %in% exomeMaleSamples]

# -------------------------------------- HapMap SNP data

library(data.table)

download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE17nnn/GSE17208/suppl/GSE17208%5FHuman660W%2DQuad%5Fv1%5F89CEU%5FFinalReport%2Etxt%2Egz", dest.file="/tmp")
download.file("ftp://webdata2:webdata2@ussd-ftp.illumina.com/MyIllumina/688e7ca8-f6af-4cde-9a33-2e990282b77f/Human660W-Quad_v1_H.csv", dest.file="/tmp")
system("gzip -d /tmp/GSE17208_Human660W-Quad_v1_89CEU_FinalReport.txt.gz")

dta <- fread("/tmp/GSE17208_Human660W-Quad_v1_89CEU_FinalReport.txt", skip=9, header=T, colClasses=c("character", "character", "character", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
setkey(dta, "Sample ID")

samples <- unique(dta$"Sample ID")

# ------------------ Split csv file into RAW Illumina FinalReports

dtb <- read.delim("/tmp/Human660W-Quad_v1_H.csv", sep=",", header=T, skip=7)
dtb.dt <- data.table(dtb)
rm(dtb)
gc()
setkey(dtb.dt, Name)
i <- 1
Chromosome <- dtb.dt[dta$"SNP Name"[(1+(657366*(i-1))):(657366*i)]]$Chr
Position <- dtb.dt[dta$"SNP Name"[(1+(657366*(i-1))):(657366*i)]]$MapInfo

dir.create("/tmp/RAW")

for (i in 1:length(samples)){
  dta.temp <- dta[(1+(657366*(i-1))):(657366*i)]
  dta.temp[, Position:=Position]
  dta.temp[, Chromosome:=Chromosome]
  write.table(dta.temp, file=file.path("tmp/RAW", paste0(dta.temp$"Sample ID"[1], "_FinalReport.txt")), quote=F, sep="\t", col.names=T, row.names=F)
}

# ------------------ Convert Illumina FinalReports to MAD files

dir.create("/tmp/rawData")

system(paste0("perl ", file.path(system.file("extdata", package = "MADloy"), "BS2MADGT.pl"), "/tmp/RAW /tmp/rawData"))