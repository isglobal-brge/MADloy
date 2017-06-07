### Download data from GEO
#
#FinalReport
#
#ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16894/suppl/GSE16894%5F1M%2DDuo%5FCEU%5FFinal%5FCall%5FReport%2Ecsv%2Egz
#gzip -d GSE16894_1M-Duo_CEU_Final_Call_Report.csv.gz
#
#Platform annot
#
#download https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE16894&format=file
#

library(data.table)

dta <- fread("GSE16894_1M-Duo_CEU_Final_Call_Report.csv", skip=, sep=",", dec=".", header=T, colClasses=c("character", "character", "character", "character", "numeric", "numeric", "numeric", "numeric","numeric","numeric","numeric"))

annot <- fread("GPL6984_Human1M-Duov3_B_noHead.csv", header=T, skip=7, stringsAsFactors=F, sep=",")
setkey(annot, Name)
dir.create("rawData")
samples <- unique(dta$"Sample ID")
for (samp in samples){
	id <- strsplit(samp, "-")[[1]][3]
	dtb <- dta[dta$"Sample ID" == samp]
	dtb <- dtb[, c(1, 10, 11), with=F]
	setnames(dtb, c("Name", "Log.R.Ratio", "B.Allele.Freq"))
	setkey(dtb, Name)
	dtb$Chr <- annot$Chr
	dtb$Position <- annot$MapInfo
	dtb$Genotype <- rep("NC", nrow(dtb))
	dtb$Genotype[ dtb$B.Allele.Freq >0.95 ] <- "BB"
	dtb$Genotype[ dtb$B.Allele.Freq <0.05 ] <- "AA"
	dtb$Genotype[ dtb$B.Allele.Freq > 0.4 &  dtb$B.Allele.Freq < 0.6 ] <- "AB"
	dtb <- dtb[, .(Name, Chr, Position, Log.R.Ratio, B.Allele.Freq, Genotype)]
	write.table(dtb, file=paste0("rawData/", id, ".txt"), quote=F, sep="\t", row.names=F, col.names=T)
}
