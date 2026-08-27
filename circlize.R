## circlize.R
## Naohiro Makise, Masahito Kawazu
## last updated 20240607

args <- commandArgs(trailingOnly = T)

## get arguments
CASE <- args[1]
SAMPLE <- args[2]
DAY <- args[3]

## merge
DIR<-paste0("/rgdata/home/makise/",CASE,"/",SAMPLE)
INDIR<-paste0(DIR,"/",DAY)
OUTDIR<-paste0(INDIR,"/circlize") 
PREF<-paste0(SAMPLE,"-",DAY,"-hacm4")

setwd(DIR)

## mkdir
if (!dir.exists(INDIR)){dir.create(INDIR)}
if (!dir.exists(OUTDIR)){dir.create(OUTDIR)}


############################################

# load library
library(dplyr)
library(tidyr)
library(toprdata)
library(circlize)
library(stringr)

#####################################


# trans3
read.table(paste0(INDIR,"/SVIM_signatures/",PREF,"-trans3.bed")) -> t3bedmm
data(ENSGENES)

# remove chrM
t3bedm<-filter(t3bedmm, !str_detect(V1, "chrM"))
t3bed<-filter(t3bedm, !str_detect(V6, "chrM"))

t3bed %>% mutate(transid = rep(c(1:(nrow(t3bed)/2)), each = 2)) %>% mutate(genename = "") -> t3bed_id
t3bed_id %>% arrange(V2) %>% arrange(V1) -> t3bed_id_ar

coln <- 1

for (i in 1:nrow(ENSGENES)){
  
  for (j in coln:nrow(t3bed_id_ar)){
    if (t3bed_id_ar[j,1] != ENSGENES[i,1]) next
    
    if (t3bed_id_ar[j,2] < ENSGENES[i,2]) next
    
    if (t3bed_id_ar[j,2] < ENSGENES[i,3]){
      t3bed_id_ar[j,8] <- ENSGENES[i,4]
      next
    }
    coln <- j
    print(coln)
    break
  }
}

t3bed_id_ar %>% arrange(transid) -> t3bed_gene

write.table(t3bed_gene, file =paste0(OUTDIR,"/",PREF,"-trans3_gene.bed"), sep ="\t")

# BND_source and BND_destination
t3beds = t3bed_gene[seq(1, nrow(t3bed_gene), by =2),]
t3beds4 <- t3beds[,1:4] %>% 
  mutate(V2 = V2 - 3000000) %>%
  mutate(V3 = V3 + 3000000)
t3bedd = t3bed_gene[seq(2, nrow(t3bed_gene), by =2),]
t3bedd4 <- t3bedd[,1:4] %>% 
  mutate(V2 = V2 - 3000000) %>%
  mutate(V3 = V3 + 3000000)

# colors
rand_color(nrow(t3bedd), transparency = 0.5) -> colorst
colorst %>% str_sub(start =1, end = -3) -> colors
colors %>% rep(each=2) -> colors2

linecolors2 <- colors2
linecolors2[t3bed_gene[,8]==""] <- "white"

# write circos-nogene
png(paste0(OUTDIR,"/",PREF,"-t3bed-circos-nogene.png"), width = 1000, height = 1000)
circos.initializeWithIdeogram(plotType = NULL)
##circos.genomicLabels(t3bed_gene, labels.column = 8, 
#                     side = "outside",
##                     cex = 1.5,
##                     col = colors2, 
##                     line_col = linecolors2)
circos.genomicIdeogram()
circos.genomicLink(t3beds4, 
                   t3bedd4, 
                   col = colorst, 
                   border = NA)
circos.clear()
dev.off()


# write circos
png(paste0(OUTDIR,"/",PREF,"-t3bed-circos.png"), width = 1000, height = 1000)
circos.initializeWithIdeogram(plotType = NULL)
circos.genomicLabels(t3bed_gene, labels.column = 8, 
                     side = "outside",
                     cex = 1.5,
                     col = colors2, 
                     line_col = linecolors2)
circos.genomicIdeogram()
circos.genomicLink(t3beds4, 
                   t3bedd4, 
                   col = colorst, 
                   border = NA)
circos.clear()
dev.off()

########################################
# trans2
read.table(paste0(INDIR,"/SVIM_signatures/",PREF,"-trans2.bed")) -> t2bedmm
data(ENSGENES)

# remove chrM
t2bedm<-filter(t2bedmm, !str_detect(V1, "chrM"))
t2bed<-filter(t2bedm, !str_detect(V6, "chrM"))



t2bed %>% mutate(transid = rep(c(1:(nrow(t2bed)/2)), each = 2)) %>% mutate(genename = "") -> t2bed_id
t2bed_id %>% arrange(V2) %>% arrange(V1) -> t2bed_id_ar

coln <- 1

for (i in 1:nrow(ENSGENES)){
  
  for (j in coln:nrow(t2bed_id_ar)){
    if (t2bed_id_ar[j,1] != ENSGENES[i,1]) next
    
    if (t2bed_id_ar[j,2] < ENSGENES[i,2]) next
    
    if (t2bed_id_ar[j,2] < ENSGENES[i,3]){
      t2bed_id_ar[j,8] <- ENSGENES[i,4]
      next
    }
    coln <- j
    print(coln)
    break
  }
}

t2bed_id_ar %>% arrange(transid) -> t2bed_gene

write.table(t2bed_gene, file =paste0(OUTDIR,"/",PREF,"-trans2_gene.bed"), sep ="\t")
