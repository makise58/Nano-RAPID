## msarc.R
## Naohiro Makise, Masahito Kawazu
## last updated 20240607

args <- commandArgs(trailingOnly = T)

## get arguments
CASE <- args[1]
SAMPLE <- args[2]
DAY <- args[3]

########################
#clean R's brain except CASE, SAMPLE, DAY
rm(list=ls()[!(ls() %in% c("CASE","SAMPLE","DAY"))])
gc()

#################################

## merge
DIR<-paste0(CASE,"/",SAMPLE)
INDIR<-paste0(DIR,"/",DAY)
OUTDIR<-paste0(INDIR,"/msarc") 
PREF<-paste0(SAMPLE,"-",DAY,"-hacm4")

setwd(DIR)

## mkdir
if (!dir.exists(INDIR)){dir.create(INDIR)}
if (!dir.exists(OUTDIR)){dir.create(OUTDIR)}

###################################
# load packages tidyverse, openxlsx
library(tidyverse)
library(openxlsx)

#####################################
# read table
target_cpg <- read.xlsx("Rdata/20230822_target_cpg_only_bed_forR.xlsx", colNames = T , rowNames = F)

# define callback
callback <- DataFrameCallback$new(function(chunk, pos) {
  chunk %>% mutate(chr_pos = paste0(X1,":",X3))->filtered
  filtered <- filtered[filtered$chr_pos %in% target_cpg$chr_pos,]
  if(pos == 0) {
    write_tsv(filtered, paste0(OUTDIR,"/",PREF,"-filtered_file.tsv"))
  } else {
    write_tsv(filtered, paste0(OUTDIR,"/",PREF,"-filtered_file.tsv"), append = TRUE, col_names = FALSE)
  }
})

# chunk size 1000000
read_tsv_chunked(paste0(DIR,"/",PREF,"-modkit.bed"), col_names = FALSE, callback, chunk_size = 1000000)


######################################################
filtered <- read.table(paste0(OUTDIR,"/",PREF,"-filtered_file.tsv"), sep = "\t", header = F)

filtered %>% rename("chr_pos" = 19) -> filtered_rev

target_cpg %>% 
  left_join(filtered_rev, by = "chr_pos") -> target_cpg_C_num

save(filtered, filtered_rev, target_cpg, target_cpg_C_num, file = paste0(OUTDIR,"/",PREF,"-filtered_target_cpg.RData"))
write.table(target_cpg_C_num, file =paste0(OUTDIR,"/",PREF,"-target_cpg_C_num.tsv"), sep ="\t", row.names = TRUE, col.names =NA)

#######################################################

target_cpg_C_num[is.na(target_cpg_C_num)] <- 0

target_cpg_C_num %>%
  group_by(cgid) %>%
  summarise(depth = sum(V10), meth = sum(V12), hmeth = sum(V14)) %>%
  mutate(beta = (meth + hmeth) / depth) -> target_cpg_C_num_summ

save(target_cpg_C_num_summ, file = paste0(OUTDIR,"/",PREF,"-target_cpg_C_num_summ.RData"))
write.table(target_cpg_C_num_summ, file =paste0(OUTDIR,"/",PREF,"-target_cpg_C_num_summ.tsv"), sep ="\t", row.names = TRUE, col.names =NA)

sarc_betas_10k_median_t <- read.table("Rdata/sarc_betas_10k_median_t.tsv", sep = "\t", header = T)

sarc_betas_10k_median_t %>% 
  left_join(target_cpg_C_num_summ, by = "cgid") %>%
  mutate(beta_med = ifelse(is.nan(beta), median, beta)) -> target_cpg_C_num_summ_med

save(target_cpg_C_num_summ_med, file = paste0(OUTDIR,"/",PREF,"-target_cpg_C_num_summ_med.RData"))
write.table(target_cpg_C_num_summ_med, file =paste0(OUTDIR,"/",PREF,"-target_cpg_C_num_summ_med.tsv"), sep ="\t", row.names = TRUE, col.names =NA)

target_cpg_C_num_summ_med %>% arrange(cgid) -> target_cpg_C_num_summ_med_ar
save(target_cpg_C_num_summ_med_ar, file = paste0(OUTDIR,"/",PREF,"-target_cpg_C_num_summ_med_ar.RData"))
write.table(target_cpg_C_num_summ_med_ar, file =paste0(OUTDIR,"/",PREF,"-target_cpg_C_num_summ_med_ar.tsv"), sep ="\t", row.names = TRUE, col.names =NA)


########################################

library(tidyverse)
library(limma)
library(openxlsx)
library(randomForest)
library(glmnet)
library(dendextend)

################################
# load datas
load("Rdata/betas_sarcref_10k.ba.RData")
load("Rdata/rf.pred.sarc.nt1000.RData")
load("Rdata/nt1000.2b.calfit.RData")

# make annotation vector
y = annosarcfinal$`Methylation.Class.Name`
y <- c(y, "UNKNOWN")

# sort ref betas
Name_list <- names(betas_sarcref_10k) %>% sort()
betas_sarcref_10k %>% select(all_of(Name_list)) -> betas_sarcref_10k_ar

# transposition
target_cpg_C_num_summ_med_ar_t <- as.matrix(t(target_cpg_C_num_summ_med_ar))


#######################################
# 1 line
s_betas <- as.matrix(target_cpg_C_num_summ_med_ar[,8])
row.names(s_betas) <- target_cpg_C_num_summ_med_ar[,1]
# rbind
betas_ref_s_10k <- rbind(betas_sarcref_10k_ar, s_betas[,1])
row.names(betas_ref_s_10k)[[1078]] <- SAMPLE

########################################
# random forest

s_betas_t <- t(s_betas)

message("running random forest ...",Sys.time())

rf.scores <- predict(rf.pred,s_betas_t[,match(rownames(rf.pred$importance),colnames(s_betas_t))],type="prob")
probs <- predict(cv.calfit$glmnet.fit,newx=rf.scores,type="response",s=cv.calfit$lambda.1se)[,,1] # use lambda estimated by 10fold CVlambda
meth_class <- as.matrix(probs) %>% as.data.frame() %>% arrange(-V1)

write.table(meth_class, file =paste0(OUTDIR,"/",PREF,"-meth_class.tsv"), sep ="\t", row.names = TRUE, col.names =NA)

save(rf.scores, probs, meth_class, file=paste0(OUTDIR,"/",PREF,"-sarc_s_RF.RData"))

message("finished random forest ...",Sys.time())
########################################
# tSNE
library(Rtsne)
library(RSpectra)
source("R/RSpectra_pca.R")

#define colors
color_map <- function(y_values) {
  hex_colors <- c('methylation class rhabdomyosarcoma (alveolar)' = "#FDA746",
                  'methylation class rhabdomyosarcoma (embryonal)' = "#F79360",
                  'methylation class synovial sarcoma' = "#7CC24B",
                  'methylation class myxoid liposarcoma' = "#F5D21A",
                  'methylation class small blue round cell tumour with BCOR alteration' = "#45565F",
                  'methylation class small blue round cell tumour with CIC alteration' = "#224571",
                  'methylation class Ewingﾂｴs sarcoma' = "#2F4B81",
                  'methylation class fibrous dysplasia' = "#5B4099",
                  'methylation class chordoma' = "#1F3761",
                  'methylation class osteosarcoma (high grade)' = "#2E3092",
                  'methylation class undifferentiated sarcoma' = "#218843",
                  'methylation class chondrosarcoma (IDH group A)' = "#6282C2",
                  'methylation class chondrosarcoma (group B)' = "#326392",
                  'methylation class chondrosarcoma (mesenchymal)' = "#4B5DAA",
                  'methylation class chondrosarcoma (group A)' = "#005DA4",
                  'methylation class extraskeletal myxoid chondrosarcoma' = "#4FB852",
                  'methylation class desmoplastic small round cell tumour' = "#44BDA4",
                  'methylation class alveolar soft part sarcoma' = "#62C7C6",
                  'methylation class sarcoma (RMS-like)' = "#F32608",
                  'methylation class solitary fibrous tumour' = "#A5257F",
                  'methylation class sarcoma (MPNST-like)' = "#B4080B",
                  'methylation class angiosarcoma' = "#CB5827",
                  'methylation class gastrointestinal stromal tumour' = "#D0A82F",
                  'methylation class leiomyosarcoma' = "#DDB081",
                  'methylation class control (reactive tissue)' = "#CACBCD",
                  'methylation class chordoma (dedifferentiated)' = "#172948",
                  'methylation class lipoma' = "#FFF450",
                  'methylation class dermatofibrosarcoma protuberans' = "#C64A9B",
                  'methylation class giant cell tumour of bone' = "#1983C3",
                  'methylation class malignant peripheral nerve sheath tumour' = "#B55A27",
                  'methylation class osteoblastoma' = "#6B66AE",
                  'methylation class desmoid-type fibromatosis' = "#DB6E9D",
                  'methylation class low-grade fibromyxoid sarcoma' = "#70003E",
                  'methylation class neurofibroma' = "#975929",
                  'methylation class neurofibroma (plexiform)' = "#6E4027",
                  'methylation class schwannoma' = "#9A5B33",
                  'methylation class infantile fibrosarcoma' = "#AB0066",
                  'methylation class rhabdomyosarcoma (MYOD1)' = "#FFB946",
                  'methylation class Kaposi sarcoma' = "#DF7D38",
                  'methylation class melanoma (cutaneous)' = "#000000",
                  'methylation class sclerosing epithelioid fibrosarcoma' = "#560030",
                  'methylation class well- / dedifferentiated liposarcoma' = "#FFC20D",
                  'methylation class control (muscle tissue)' = "#7C7E82",
                  'methylation class malignant rhabdoid tumour' = "#55A458",
                  'methylation class chondrosarcoma (clear cell)' = "#004892",
                  'methylation class chondroblastoma' = "#A3B1DA",
                  'methylation class chondrosarcoma (IDH group B)' = "#1B75BC",
                  'methylation class squamous cell carcinoma (cutaneous)' = "#231F20",
                  'methylation class endometrial stromal sarcoma (low grade)' = "#65356C",
                  'methylation class atypical fibroxanthoma / pleomorphic dermal sarcoma' = "#B6D995",
                  'methylation class clear cell sarcoma of the kidney' = "#C98DBE",
                  'methylation class epithelioid sarcoma' = "#6ABE4F",
                  'methylation class myositis ossificans' = "#E8A9CC",
                  'methylation class epithelioid haemangioendothelioma' = "#FE8C17",
                  'methylation class inflammatory myofibroblastic tumour' = "#901E5B",
                  'methylation class nodular fasciitis' = "#FACAD1",
                  'methylation class angiomatoid fibrous histiocytoma' = "#9CD089",
                  'methylation class myositis proliferans' = "#C57D95",
                  'methylation class ossifying fibromyxoid tumour' = "#82C66F",
                  'methylation class Langerhans cell histiocytosis' = "#4D2D80",
                  'methylation class leiomyoma' = "#DED49D",
                  'methylation class clear cell sarcoma of soft parts' = "#49C1BA",
                  'methylation class angioleiomyoma / myopericytoma' = "#F36E2B",
                  'methylation class endometrial stromal sarcoma (high grade)' = "#3C1448",
                  'methylation class control (blood)' = "#67686C",
                  'UNKNOWN' = "#111111"
  )
  return(sapply(y_values, function(val) hex_colors[as.character(val)]))
}

message("running tSNE ...",Sys.time())

# calculate first 65 PCs
pca <- prcomp_svds(betas_ref_s_10k,k=65)

# calculate tSNE
res15 <- Rtsne(pca$x,pca=F,perplexity =15,max_iter=3000,theta=0,verbose=T)
res30 <- Rtsne(pca$x,pca=F,perplexity =30,max_iter=3000,theta=0,verbose=T)

# scatterplot tSNE
tsne_path <- paste0(OUTDIR,"/",PREF,"-tsne15.png")
png(tsne_path, width = 1000, height = 1000)
plot(res15$Y,pch=19,col=color_map(y))
dev.off()

tsne_path <- paste0(OUTDIR,"/",PREF,"-tsne30.png")
png(tsne_path, width = 1000, height = 1000)
plot(res30$Y,pch=19,col=color_map(y))
dev.off()

# calculate distance
res15_mat <- res15[["Y"]]
sx <- res15_mat[1078,1]
sy <- res15_mat[1078,2]
res15_mat %>% 
  cbind(y) %>%
  as.data.frame %>%
  mutate(dist = sqrt((as.numeric(V1)-sx)^2+(as.numeric(V2)-sy)^2)) %>%
  arrange(dist) -> res15_dist
write.table(res15_dist, file =paste0(OUTDIR,"/",PREF,"-res15_dist.tsv"), sep ="\t", row.names = TRUE, col.names =NA)

res30_mat <- res30[["Y"]]
sx <- res30_mat[1078,1]
sy <- res30_mat[1078,2]
res30_mat %>% 
  cbind(y) %>%
  as.data.frame %>%
  mutate(dist = sqrt((as.numeric(V1)-sx)^2+(as.numeric(V2)-sy)^2)) %>%
  arrange(dist) -> res30_dist
write.table(res30_dist, file =paste0(OUTDIR,"/",PREF,"-res30_dist.tsv"), sep ="\t", row.names = TRUE, col.names =NA)

save(res15, res15_dist, res30, res30_dist, pca, file=paste0(OUTDIR,"/",PREF,"-tsne.RData"))

message("finished tSNE ...",Sys.time())

##############################################
#HC, answer&color, answer&black

message("start hierarchical clustering ...",Sys.time())

# Eucledian distance
message("calculating Eucledian distance ...",Sys.time())
dist_mat <- dist(as.data.frame(betas_ref_s_10k))

save(dist_mat, file=paste0(OUTDIR,"/",PREF,"-Eucledian_dist.RData"))


###############################################
# method = ward.d2
hclust_ward <- hclust(dist_mat, method = "ward.D2")
dend_ward <- as.dendrogram(hclust_ward)

# ggdend
gdend_ward <- as.ggdend(dend_ward)
gdend_ward$labels$angle <- 90
#gdend_ward$labels <- gdend_ward$labels %>% mutate(col = color_map(y))

# add colors
gdend_ward$labels %>% left_join(sarc_GSM_to_col_final, by ="label") %>%
  select(1,2,5,6,7,8) %>%
  rename("label" = 5, "col" = 6) ->gdend_ward$labels

# draw png
pngpath <- paste0(OUTDIR,"/",PREF,"-sarc_hclust_ward.png")
png(file = pngpath, width = 19000, height = 10000)
ggplot(gdend_ward)
dev.off()

message("finished hierarchical clustering ...",Sys.time())

