## QDNAseq.R
## Naohiro Makise, Masahito Kawazu 
## last updated 20240422

args <- commandArgs(trailingOnly = T)

## get arguments
CASE <- args[1]
SAMPLE <- args[2]
DAY <- args[3]

## merge
DIR<-paste0("/rgdata/home/makise/",CASE,"/",SAMPLE)
INDIR<-paste0(DIR,"/",DAY)
OUTDIR<-paste0(INDIR,"/QDNAseq") 
PREF<-paste0(SAMPLE,"-",DAY,"-hacm4")

setwd(DIR)

## mkdir
if (!dir.exists(INDIR)){dir.create(INDIR)}
if (!dir.exists(OUTDIR)){dir.create(OUTDIR)}


#########################################
# load library
library(QDNAseq)
library(QDNAseq.hg38)

# set BAM
BAM<-paste0(PREF,"-u1000.bam")

# bin: 1,5,10,15,30,50,100,500,1000
binvec <- c(30,50,100,500,1000)

for (i in 1:5){
  binset <- binvec[i]
  
  binset
  file_name <- paste(PREF, "-u1000_bin", sprintf("%04d", binset), sep ="")
  
  
  bins<- getBinAnnotations(binSize = binset, genome ="hg38", type = "SR50")
  
  readCounts <- binReadCounts(bins, bamfiles = BAM)
  
  ## plot
  plot(readCounts, logTransform=FALSE)
  
  readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
  ## GC content vs median read counts縲
  isobarPlot(readCountsFiltered)
  
  readCountsFiltered <- estimateCorrection(readCountsFiltered)
  ## noiseplot
  noisePlot(readCountsFiltered)
  
  ## normalize, smoothen
  copyNumbers <- correctBins(readCountsFiltered)
  copyNumbersNormalized <- normalizeBins(copyNumbers)
  copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
  plot(copyNumbersSmooth)
  
  ##export
  exportBins(copyNumbersSmooth, file= paste(OUTDIR, "/", file_name, ".txt", sep =""))
  exportBins(copyNumbersSmooth, file= paste(OUTDIR, "/", file_name, ".igv", sep =""), format = "igv")
  exportBins(copyNumbersSmooth, file= paste(OUTDIR, "/", file_name, ".bed", sep =""), format = "bed")
  
  #segmentation
  copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt")
  copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
  
  png(paste(OUTDIR, "/", file_name, "_seg.png", sep= ""))
  plot(copyNumbersSegmented)
  dev.off()
  
  pdf(paste(OUTDIR, "/", file_name, "_seg.pdf", sep= ""))
  plot(copyNumbersSegmented)
  dev.off()
 
}
