#Naohiro Makise, Masahito Kawazu

#sarcoma methyl classify
max.print = 1000
#options(stringsAsFactors = FALSE)
#options(scipen = 999)
rm(list=ls())


##preprocessing


library(minfi)
library(GEOquery)
library(limma)
library(openxlsx)
library(stringr)

source(file.path("R","MNPprocessIDAT_functions.R"))

dir.create("results")

## get sample annotation from GEO
gse <- getGEO("GSE140686", GSEMatrix=TRUE, getGPL=FALSE)
annoEPIC <- pData(gse$`GSE140686-GPL21145_series_matrix.txt.gz`)
write.table(annoEPIC, file="GSE140686-GPL21145_series_matrix.tsv", sep = "\t", row.names = T, col.names = NA)
annoEPICref <- subset(annoEPIC, str_detect(annoEPIC$title,"reference"))
write.table(annoEPICref, file="GSE140686-GPL21145_series_ref_matrix.tsv", sep = "\t", row.names = T, col.names = NA)
annoEPICval <- subset(annoEPIC, str_detect(annoEPIC$title,"validation"))
write.table(annoEPICval, file="GSE140686-GPL21145_series_val_matrix.tsv", sep = "\t", row.names = T, col.names = NA)
anno450k <-  pData(gse$`GSE140686-GPL13534_series_matrix.txt.gz`)
write.table(anno450k, file="GSE140686-GPL13534_series_matrix.tsv", sep = "\t", row.names = T, col.names = NA)


annosarcref <- rbind(anno450k, annoEPICref)
write.table(annosarcref, file="GSE140686-GPL13534-GPL21145_ref_series_matrix.tsv", sep = "\t", row.names = T, col.names = NA)

#######################################################################
# 450k first
# all cases are references
# read raw data downloaded from GEO and extracted in GSE140686_RAW
filepath <- file.path("GSE140686_RAW_idat",gsub("_Grn.*","",gsub(".*suppl/","",anno450k$supplementary_file)))
RGset450k <- read.metharray(filepath,verbose=TRUE)

# Illumina normalization
message("running normalization ...",Sys.time())
Mset450k <- MNPpreprocessIllumina(RGset450k)

# probe filtering
message("probe filtering ...",Sys.time())
amb.filter <- read.table(file.path("filter","amb_3965probes.vh20151030.txt"),header=F)
epic.filter <- read.table(file.path("filter","epicV1B2_32260probes.vh20160325.txt"),header=F)
snp.filter <- read.table(file.path("filter","snp_7998probes.vh20151030.txt"),header=F)
xy.filter <- read.table(file.path("filter","xy_11551probes.vh20151030.txt"),header=F)
rs.filter <- grep("rs",rownames(Mset450k))
ch.filter <- grep("ch",rownames(Mset450k))

# filter CpG probes
remove <- unique(c(match(amb.filter[,1], rownames(Mset450k)),
                   match(epic.filter[,1], rownames(Mset450k)),
                   match(snp.filter[,1], rownames(Mset450k)),
                   match(xy.filter[,1], rownames(Mset450k)),
                   rs.filter,
                   ch.filter))

Mset450k_filtered <- Mset450k[-remove,]

# 428799 CpGs after filtering
Mset450k_filtered@NAMES
write.table(Mset450k_filtered@NAMES, file="CpG_filtered_428799.tsv", sep = "\t", col.names = F, row.names = F)


save(Mset450k_filtered,anno450k,annoEPICref,annoEPICval,file=file.path("results","Mset450k_filtered.RData"))  

rm(Mset)
gc()


###################################################
# Next EPIC cohort only
# read raw data downloaded from GEO and extracted in GSE140686_RAW
filepath <- file.path("GSE140686_RAW_idat",gsub("_Grn.*","",gsub(".*suppl/","",annoEPICref$supplementary_file)))
RGsetEPICref <- read.metharray(filepath,verbose=TRUE,force=TRUE)

# Illumina normalization
message("running normalization ...",Sys.time())
MsetEPICref <- MNPpreprocessIllumina(RGsetEPICref)

save(MsetEPICref,anno450k,annoEPICref,annoEPICval,file=file.path("results","MsetEPICref.RData"))  


####################################################
# combine 450k and EPIC
Mset_sarcref <- combineArrays(Mset450k_filtered, MsetEPICref, verbose = TRUE)
# done!428230 probes!
save(Mset_sarcref,anno450k,annoEPICref,annoEPICval,annosarcref,file=file.path("results","Mset_sarcref.RData"))

Mset_sarcref@NAMES
# write 428230 probes
# write.table(Mset_sarcref@NAMES, file="CpG_filtered_428230.tsv", sep = "\t", col.names = F, row.names = F)
## done!

rm(Mset450k_filtered)
rm(MsetEPICref)
gc()

####################################################
# batch FFPE vs KRYO and batch 450k vs EPIC

message("performing batchadjustment ...",Sys.time())

methy <- getMeth(Mset_sarcref)
unmethy <- getUnmeth(Mset_sarcref)
rm(Mset_sarcref)
gc()

# get FFPE/Frozen type
ffpe <- annosarcref$`material preparation:ch1`
batch <- ifelse(ffpe == "FFPE", 2, 1)

# get 450k/EPIC type
EPIC <- annosarcref$`platform_id`
batch2 <- ifelse(EPIC == "GPL21145", 2, 1)

# remove batch effects by linear model
methy.ba <- 2^removeBatchEffect(log2(methy +1), batch, batch2)
unmethy.ba <- 2^removeBatchEffect(log2(unmethy +1), batch, batch2)

# extract effects to adjust diagnostic samples
s.kryo.450k <- min(which(batch == 1 & batch2 == 1))
s.ffpe.450k <- min(which(batch == 2 & batch2 == 1))
s.kryo.EPIC <- min(which(batch == 1 & batch2 == 2))
s.ffpe.EPIC <- min(which(batch == 2 & batch2 == 2))

methy.coef <- unmethy.coef <- list()
methy.coef[["KRYO.450K"]] <- log2(methy.ba[, s.kryo.450k]) - log2(methy[, s.kryo.450k] +1)
methy.coef[["FFPE.450K"]] <- log2(methy.ba[, s.ffpe.450k]) - log2(methy[, s.ffpe.450k] +1)
methy.coef[["KRYO.EPIC"]] <- log2(methy.ba[, s.kryo.EPIC]) - log2(methy[, s.kryo.EPIC] +1)
methy.coef[["FFPE.EPIC"]] <- log2(methy.ba[, s.ffpe.EPIC]) - log2(methy[, s.ffpe.EPIC] +1)

unmethy.coef[["KRYO.450K"]] <- log2(unmethy.ba[, s.kryo.450k]) - log2(unmethy[, s.kryo.450k] +1)
unmethy.coef[["FFPE.450K"]] <- log2(unmethy.ba[, s.ffpe.450k]) - log2(unmethy[, s.ffpe.450k] +1)
unmethy.coef[["KRYO.EPIC"]] <- log2(unmethy.ba[, s.kryo.EPIC]) - log2(unmethy[, s.kryo.EPIC] +1)
unmethy.coef[["FFPE.EPIC"]] <- log2(unmethy.ba[, s.ffpe.EPIC]) - log2(unmethy[, s.ffpe.EPIC] +1)

## write coef
write.table(methy.coef[["KRYO.450K"]], sep = "\t", row.names = T, col.names = F, file = "sarc_kryo_450k_met_coef.tsv")
write.table(methy.coef[["FFPE.450K"]], sep = "\t", row.names = T, col.names = F, file = "sarc_ffpe_450k_met_coef.tsv")
write.table(methy.coef[["KRYO.EPIC"]], sep = "\t", row.names = T, col.names = F, file = "sarc_kryo_epic_met_coef.tsv")
write.table(methy.coef[["FFPE.EPIC"]], sep = "\t", row.names = T, col.names = F, file = "sarc_ffpe_epic_met_coef.tsv")

write.table(unmethy.coef[["KRYO.450K"]], sep = "\t", row.names = T, col.names = F, file = "sarc_kryo_450k_un_coef.tsv")
write.table(unmethy.coef[["FFPE.450K"]], sep = "\t", row.names = T, col.names = F, file = "sarc_ffpe_450k_un_coef.tsv")
write.table(unmethy.coef[["KRYO.EPIC"]], sep = "\t", row.names = T, col.names = F, file = "sarc_kryo_epic_un_coef.tsv")
write.table(unmethy.coef[["FFPE.EPIC"]], sep = "\t", row.names = T, col.names = F, file = "sarc_ffpe_epic_un_coef.tsv")


# save batch effects 
save(methy.coef,unmethy.coef,file=file.path("results","sarcref.ba.coef.RData"))

# recalculate betas, illumina like
betas_sarcref <- methy.ba / (methy.ba +unmethy.ba +100)
betas_sarcref <- as.data.frame(t(betas_sarcref))
save(betas_sarcref,anno450k,annoEPICref,annoEPICval,annosarcref,file=file.path("results","betas_sarcref.ba.RData"))  
message("preprocessing finished ...",Sys.time())



#################################################

##bind annosarcref and methylation class xlsx
library(dplyr)
annosarcclass <- read.xlsx("41467_2020_20603_MOESM4_ESM.xlsx", colNames = TRUE, rowNames = FALSE, check.names = FALSE)
annosarcfinal <- left_join(annosarcref, annosarcclass, by=c(`description.1`="ID"))
 
write.table(annosarcfinal, file="GSE140686-GPL13534-GPL21145_ref_class_series_matrix.tsv", sep = "\t", row.names = T, col.names = NA)


##tSNE
rm(list=ls())

library(Rtsne)
library(RSpectra)

source(file.path("R","RSpectra_pca.R"))

message("loading preprocessed data ...",Sys.time())
load(file.path("results","betas_sarcref.ba.RData"))

# methylation classes
annosarcfinal <- read.table("GSE140686-GPL13534-GPL21145_ref_class_series_matrix.tsv", check.names = F, header = T, row.names =1)
y <- as.factor(annosarcfinal$`Methylation.Class.Name`)

#####################################################################
# sd filtering to 10k probes
betas_sarcref_10k <- betas_sarcref_15k[,order(-apply(betas_sarcref_15k,2,sd))[1:10000]]
rm(betas_sarcref_15k)
gc()
# write 10k probes
write.table (t(betas_sarcref_10k[1,1:10000]), sep = "\t", row.names = T, col.names = NA, file = "sarcref_beta10000.tsv")

# calculate first 65 PCs
pca <- prcomp_svds(betas_sarcref_10k,k=65)

# calculate tSNE
res15 <- Rtsne(pca$x,pca=F,perplexity=15,max_iter=3000,theta=0,verbose=T)
scatterplot tSNE
plot(res15$Y,pch=19,col=y, main = "perplexity15, iteration3000")

# calculate tSNE
res30 <- Rtsne(pca$x,pca=F,perplexity=30,max_iter=3000,theta=0,verbose=T)
scatterplot tSNE
plot(res30$Y,pch=19,col=y, main = "perplexity30, iteration3000")

save(betas_sarcref_10k,annosarcfinal,res15,res30,file=file.path("results","betas_sarcref_10k.ba.RData"))

#####################################################################


############################################
# cross validation sarcoma
max.print = 1000
#options(stringsAsFactors = FALSE)
#options(scipen = 999)
rm(list=ls())

library(randomForest)
library(parallel)
library(minfi)
library(limma)

ntrees <- 1000 # 10000 in paper
cores <- 4
seed <- 180314
p <- 10000
folds <- 3

message("loading filtered Mset ...",Sys.time())
load(file.path("results","Mset_sarcref.RData"))

annosarcfinal <- read.table("GSE140686-GPL13534-GPL21145_ref_class_series_matrix.tsv", check.names = F, header = T, row.names =1)
y <- as.factor(annosarcfinal$`Methylation.Class.Name`)
batch1 <- as.factor(annosarcfinal$`material preparation:ch1`) # FFPE or KRYO
batch2 <- as.factor(annosarcfinal$`platform_id`) # GPL21145(EPIC) or GPL13534(450k)
batch12 <- as.factor(paste0(batch1, ".", batch2))

source(file.path("R","makefolds.R"))
source(file.path("R","train.R"))
source(file.path("R","calculateCVfold_sarc_2batch.R"))
source(file.path("R","batchadjust_sarc_2batch.R"))

if(!file.exists(file.path("CV","nfolds.RData"))){
  dir.create("CV",showWarnings = FALSE)
  nfolds <- makenestedfolds(y,folds)
  save(nfolds,file=file.path("CV","nfolds.RData"))
}
load(file.path("CV","nfolds.RData"))

message("performing nested CV ...", Sys.time())
message("check minimal class sizes for inner training loops")

# check minimal class sizes for inner training loops
minclasssize <- matrix(0,ncol=length(nfolds),nrow=length(nfolds))
for(i in 1:length(nfolds)){
  for(j in 1:length(nfolds))
    minclasssize[i,j]  <- min(table(y[nfolds[[i]][[2]][[j]]$train]))
}
colnames(minclasssize) <- paste0("innfold",1:folds)
rownames(minclasssize) <- paste0("fold",1:folds)
print(minclasssize)

for(K in 1:folds){
  
  for(k in 0:folds){
    
    if(k>0){  message("calculateing fold ",K,".",k,"  ...",Sys.time())
      fold <- nfolds[[K]][[2]][[k]]
    }else{
      message("calculateing outer fold ",K,"  ...",Sys.time())
      fold <- nfolds[[K]][[1]][[1]]
    }
    
    rf.scores <- calcultateCVfold_sarc_2batch(Mset_sarcref,y,batch1,batch2,batch12,fold,p,cores,ntrees)
    
    fname <- paste("CVfold",K,k,"RData",sep=".")
    save(rf.scores,file=file.path("CV",fname))
    
    rm(rf.scores)
    gc()
  }
}
message("finished ...",Sys.time())

#################################
# Calibration
# nt1000.2b
library(rmarkdown)
library(glmnet)
library(doParallel)
library(HandTill2001)

cores <- 4 # default 4

registerDoParallel(cores)

message("loading data ...",Sys.time())
load(file.path("results","betas_sarcref_10k.ba.RData")) 
y <- as.factor(annosarcfinal$`Methylation.Class.Name`)

load(file.path("nt1000.2b.CV","nt1000.2b.nfolds.RData"))

for(i in 1:length(nfolds)){
  scores <- list() 
  idx <- list()
  for(j in 1:length(nfolds)){
    fname <- paste0("nt1000.2b.CVfold.",i,".",j,".RData")
    load(file.path("nt1000.2b.CV",fname))
    scores[[j]] <- rf.scores
    idx[[j]] <- nfolds[[i]][[2]][[j]]$test
  }
  scores <- do.call(rbind,scores)
  idx <- unlist(idx)
  y <- annosarcfinal$`Methylation.Class.Name`[idx]         
  
  message("fitting calbriation model fold ",i," ...",Sys.time())
  # fit multinomial logistic ridge regression model
  suppressWarnings(cv.calfit <- cv.glmnet(y=y,x=scores,family="multinomial",type.measure="mse",
                                          alpha=0,nlambda=100,lambda.min.ratio=10^-6,parallel=TRUE))
  
  fname <- paste0("nt1000.2b.CVfold.",i,".",0,".RData")
  load(file.path("nt1000.2b.CV",fname))
  
  message("calibrating raw scores fold ",i," ...",Sys.time())
  probs <- predict(cv.calfit$glmnet.fit,newx=rf.scores,type="response"
                   ,s=cv.calfit$lambda.1se)[,,1] # use lambda estimated by 10fold CVlambda
  
  
  err <- sum(colnames(probs)[apply(probs,1,which.max)] != annosarcfinal$`Methylation.Class.Name`[nfolds[[i]][[1]][[1]]$test])/length(nfolds[[i]][[1]][[1]]$test)
  
  message("misclassification error: ",err)
  
  fname_probs <- paste0("nt1000.2b.probsCVfold.",i,".",0,".RData")
  save(probs,file=file.path("nt1000.2b.CV",fname_probs))
}

scores <- list()
idx <- list()
for(i in 1:length(nfolds)){
  fname <- paste0("nt1000.2b.CVfold.",i,".",0,".RData")
  load(file.path("nt1000.2b.CV",fname))
  scores[[i]] <- rf.scores
  idx[[i]] <- nfolds[[i]][[1]][[1]]$test
}
scores <- do.call(rbind,scores)

probl <- list()
for(i in 1:length(nfolds)){
  fname <- paste0("nt1000.2b.probsCVfold.",i,".",0,".RData")
  load(file.path("nt1000.2b.CV",fname))
  probl[[i]] <- probs
}
probs <- do.call(rbind,probl)


idx <- unlist(idx)
y <- annosarcfinal$`Methylation.Class.Name`[idx] 

ys <- colnames(scores)[apply(scores,1,which.max)]
yp <- colnames(probs)[apply(probs,1,which.max)]

errs <- sum(y!=ys)/length(y)
errp <- sum(y!=yp)/length(y)

message("overall misclassification error scores: ",errs)
message("overall misclassification error calibrated: ",errp)

message("fitting final calibration model ...",Sys.time())

suppressWarnings(cv.calfit <- cv.glmnet(y=y,x=scores,family="multinomial",type.measure="mse",
                                        alpha=0,nlambda=100,lambda.min.ratio=10^-6,parallel=TRUE))

save(cv.calfit,file=file.path("results","nt1000.2b.calfit.RData"))

save(scores,probs,y,ys,yp,file=file.path("results","nt1000.2b.CVresults.RData"))

message("generating report ...",Sys.time())
rmarkdown::render("nt1000.2b.CVresults-fixed.Rmd")
message("finished ...",Sys.time())

###########################################################
rm(list=ls())
gc()

library(tidyverse)

# sarc_beta_median
sarc_betas_10k_m <- betas_sarcref_10k[1,]

rownames(sarc_betas_10k_m)[1] <- "median"
sarc_betas_10k_median <- sarc_betas_10k_m

for (i in 1:10000){
  if (i / 100 - round(i/100) == 0){ print(i)}
  sarc_betas_10k_median[i] <- median(betas_sarcref_10k[,i])
}

save(sarc_betas_10k_median,file=file.path("results","sarc_betas_10k_median.RData"))
write.table(sarc_betas_10k_median, file="sarc_betas_10k_median.tsv", sep = "\t", row.names = T, col.names = NA)

######################################################
#validate on validation set of 428 sarcomas

rm(list=ls())

library(minfi)
library(GEOquery)
library(limma)
library(openxlsx)

source(file.path("R","MNPprocessIDAT_functions.R"))


# get sample annotation from GEO
anno_sarc_val <- read.table("tables/GSE140686-GPL21145_series_val_matrix.tsv", check.names = F, header = T)
# yv <- anno_sarc_val$`Methylation.Class.Name`


# read raw data downloaded from GEO and extracted in GSE90496_RAW
filepath <- file.path("GSE140686_RAW_idat",gsub("_Grn.*","",gsub(".*suppl/","",anno_sarc_val$supplementary_file)))
RGset <- read.metharray(filepath,force=TRUE,verbose=TRUE)
save(RGset,anno_sarc_val,file=file.path("results","RGset.RData"))  

# Illumina normalization
message("running normalization ...",Sys.time())
Mset <- MNPpreprocessIllumina(RGset)
save(Mset,anno_sarc_val,file=file.path("results","Mset.RData"))  
rm(RGset)

# probe filtering
cpg_filtered_428230 <- read.table("tables/CpG_filtered_428230.tsv", check.names = F)

message("probe filtering ...",Sys.time())
Mset_filtered <- Mset[rownames(Mset) %in% cpg_filtered_428230[,1], ]

save(Mset_filtered,anno_sarc_val,file=file.path("results","Mset_filtered.RData"))  

rm(Mset)
gc()

#batch adjustment
message("performing batchadjustment ...",Sys.time())

methy <- getMeth(Mset_filtered)
unmethy <- getUnmeth(Mset_filtered)
rm(Mset_filtered)
gc()

# get FFPE/Frozen type
batch <- paste0(anno_sarc_val$`material preparation:ch1`,".EPIC")
# gsub("KRYO", "Frozen", batch)

# load coef
load(file.path("datas","sarcref.ba.coef.RData"))

# adjustment
methy.b <- log2(methy +1) + matrix(unlist(methy.coef[match(batch,names(methy.coef))]),ncol=length(batch))
unmethy.b <- log2(unmethy +1) + matrix(unlist(unmethy.coef[match(batch,names(unmethy.coef))]),ncol=length(batch))
methy.b[methy.b < 0] <- 0
unmethy.b[unmethy.b < 0] <- 0
methy.ba <- 2^methy.b
unmethy.ba <- 2^unmethy.b
# illumina-like beta values
betas.ba <- methy.ba / (methy.ba +unmethy.ba +100)
betas.ba <- as.data.frame(t(betas.ba))


save(betas.ba,anno_sarc_val,batch,cpg_filtered_428230,file=file.path("results","sarc_val_betas.ba.RData"))

# 10k probes
sarcref_10k <- read.table("tables/sarcref_beta10000.tsv", check.names = F)
betas.ba_10k <- betas.ba[colnames(betas.ba) %in% sarcref_10k_matr[,1]]

save(betas.ba_10k,anno_sarc_val,batch,cpg_filtered_428230,sarcref_10k,file=file.path("results","sarc_val_10k_betas.ba.RData"))
rm(betas.ba)


########################################

# random forest
load(file.path("datas","rf.pred.sarc.nt1000.RData"))
load(file.path("datas","nt1000.2b.calfit.RData"))

message("running random forest ...",Sys.time())

rf.scores <- predict(rf.pred,betas.ba_10k[,match(rownames(rf.pred$importance),colnames(betas.ba_10k))],type="prob")
probs <- predict(cv.calfit$glmnet.fit,newx=rf.scores,type="response",s=cv.calfit$lambda.1se)[,,1] # use lambda estimated by 10fold CVlambda
meth_class <- colnames(probs)[apply(probs,1,which.max)]

save(rf.scores, probs, meth_class,file=file.path("results","sarc_val_RF.RData"))


message("finished random forest ...",Sys.time())

# validation
sarc_val_paper <- read.xlsx("tables/41467_2020_20603_MOESM6_ESM.xlsx", check.names = F, colNames = T, rowNames = T)
sarc_val_anno <- read.xlsx("tables/41467_2020_20603_MOESM5_ESM_mod.xlsx", check.names = F, colNames = T, rowNames = F)
sarc_val_anno %>% rename("meth_class" = 1) %>% select(1:2) -> sarc_val_anno

sarc_val_paper %>% cbind(meth_class) %>%
                    left_join(sarc_val_anno, by = "meth_class" ) %>%
                    cbind(rowMax(probs)) %>%
                    rename("cal_score" = 18) ->sarc_val_kotae

sarc_val_kotae %>% mutate(unc_or_not = ((V12.2_MaxCalScore - 0.9) * (cal_score-0.9))>0) %>%
                  mutate(class_or_not = (Methylation.Class.Name.Abbreviated == V12.2_MaxCalDiag)) ->sarc_val_kotaeawase

table(sarc_val_kotaeawase$unc_or_not)
table(sarc_val_kotaeawase$class_or_not)

save(sarc_val_kotaeawase,file=file.path("results","sarc_val_kotaeawase.RData"))
write.table(sarc_val_kotaeawase, file="tables/sarc_val_kotaeawase.tsv", sep = "\t", row.names = T, col.names = NA)


#########################################
rm(list=ls())
gc()

# bind matrix and tSNE, HC
load(file.path("datas","betas_sarcref_10k.ba.RData"))
load(file.path("results","sarc_val_10k_betas.ba.RData"))

betas_bind_10k <- rbind(betas_sarcref_10k, betas.ba_10k)

save(betas_bind_10k,file=file.path("results","betas_bind_10k.RData"))
write.table(betas_bind_10k, file="tables/sarc_ref_val_betas_bind_10k.tsv", sep = "\t", row.names = T, col.names = NA)

#########################################

rm(res, res15, res30)
load(file.path("results","sarc_val_kotaeawase.RData"))
#tSNE, answer&color, answer&black

library(Rtsne)
library(RSpectra)

source(file.path("R","RSpectra_pca.R"))

yr <- as.factor(annosarcfinal$`Methylation.Class.Name`)
yv <- as.factor(sarc_val_kotaeawase$`meth_class`)
y <- c(yr,yv)

rvr <- rep("ref", 1077)
rvv <- rep("val", 428)
rv <- c(rvr, rvv)

#濶ｲ繧呈欠螳・color_map <- function(y_values) {
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

# calculate first 94 PCs
pca <- prcomp_svds(betas_bind_10k,k=65)

# calculate tSNE
res15 <- Rtsne(pca$x,pca=F,perplexity = 15, max_iter=3000,theta=0,verbose=T)
res30 <- Rtsne(pca$x,pca=F,perplexity = 30, max_iter=3000,theta=0,verbose=T)

# scatterplot tSNE
# first color
png("results/sarc_ref_val_p15_color.png", width = 1900, height = 1000)
plot(res15$Y,pch=19,col=color_map(y))
dev.off()

pdf("results/sarc_ref_val_p15_color.pdf", width = 40, height = 20)
plot(res15$Y,pch=19,col=color_map(y))
dev.off()

png("results/sarc_ref_val_p30_color.png", width = 1900, height = 1000)
plot(res30$Y,pch=19,col=color_map(y))
dev.off()

pdf("results/sarc_ref_val_p30_color.pdf", width = 40, height = 20)
plot(res30$Y,pch=19,col=color_map(y))
dev.off()

# secound reference black
png("results/sarc_ref_val_p15_black.png", width = 1900, height = 1000)
plot(res15$Y,pch=19,col= ifelse(rv =="ref", "black", color_map(y)))
dev.off()

pdf("results/sarc_ref_val_p15_black.pdf", width = 40, height = 20)
plot(res15$Y,pch=19,col= ifelse(rv =="ref", "black", color_map(y)))
dev.off()


png("results/sarc_ref_val_p30_black.png", width = 1900, height = 1000)
plot(res30$Y,pch=19,col= ifelse(rv =="ref", "black", color_map(y)))
dev.off()

pdf("results/sarc_ref_val_p30_black.pdf", width = 40, height = 20)
plot(res30$Y,pch=19,col= ifelse(rv =="ref", "black", color_map(y)))
dev.off()



save(betas_bind_10k,y,rv,pca,res15,res30,color_map,file=file.path("results","sarc_ref_val_tSNE.RData"))

#########################################

#HC, answer&color, answer&black

# Eucledian distance
message("calcurting Eucledian distance ...",Sys.time())
dist_mat <- dist(as.data.frame(betas_bind_10k))
save(dist_mat, betas_bind_10k, file = file.path("results", "sarc_ref_val_Eucl_dist.RData"))

##########################################
# method = complete
hclust_comp <- hclust(dist_mat, method = "complete")
dend_comp <- as.dendrogram(hclust_comp)

# ggdend
gdend_comp <- as.ggdend(dend_comp)
gdend_comp$labels$angle <- 90
#gdend_comp$labels <- gdend_comp$labels %>% mutate(col = colo_map(y))


# make table GSM to color
annosarcfinal_2 <- cbind(annosarcfinal[,2], annosarcfinal[,42])
anno_sarc_val_2 <- cbind(anno_sarc_val[,3], sarc_val_kotaeawase[,16])
anno_ref_val <- rbind(annosarcfinal_2, anno_sarc_val_2)
colnames(anno_ref_val) <- c("geo_accession", "Methylation.Class.Name")

betas_bind_10k %>% cbind(anno_ref_val) %>%
  cbind(color_map(y)) %>%
  select(10001:10003) -> sarc_GSM_to_col

sarc_GSM_to_col %>% cbind(substr(row.names(sarc_GSM_to_col),1,10)) %>% cbind(row.names(sarc_GSM_to_col)) -> sarc_GSM_to_col_2
sarc_GSM_to_col_2 %>% rename("geo_accession" =1, "hex_col"=3, "label2" = 4, "label" =5) ->sarc_GSM_to_col_3
sarc_GSM_to_col_3 %>% mutate(equal = ifelse(geo_accession == label2, TRUE, FALSE)) ->sarc_GSM_to_col_4

# confirm equal
table(sarc_GSM_to_col_4$equal)
sarc_GSM_to_col_4 %>% select(2,3,5) ->sarc_GSM_to_col_final

# add colors
gdend_comp$labels %>% left_join(sarc_GSM_to_col_final, by ="label") %>%
  select(1,2,5,6,7,8) %>%
  rename("label" = 5, "col" = 6) ->gdend_comp$labels


save(sarc_GSM_to_col_final, hclust_comp, dend_comp, gdend_comp, file = file.path("results", "sarc_ref_val_HC_comp.RData"))

# draw png
pngpath <- "results/sarc_hclust_comp.png"
png(file = pngpath, width = 19000, height = 10000)
ggplot(gdend_comp)
dev.off()


#########################################
# method = complete and circlize
gdend_comp_circle <- gdend_comp
gdend_comp_circle$labels <- gdend_comp$labels %>% mutate(angle = (90 - 360 * x /(max(x) + 1)))

pngpath <- "results/sarc_hclust_comp_circ.png"
png(file = pngpath, width = 20000, height = 20000)
ggplot() +
  geom_segment(data = gdend_comp$segments, aes(x =x, y = y, xend = xend, yend = yend, col = col)) +
  geom_text(data = gdend_comp_circle$labels, aes(x =x, y = y-10, label = label, col = col, angle = angle),hjust = 0) +
  scale_x_continuous(expand = c(0,1))+
  scale_y_reverse(expand = c(0,100)) +
  coord_polar(theta="x") +
  theme_minimal()+ 
  theme(legend.position = "none",
        axis.line.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.title.y=element_blank()) + 
  scale_color_identity()
dev.off()


#########################################


###############################################
# method = ward.d2
hclust_ward <- hclust(dist_mat, method = "ward.D2")
dend_ward <- as.dendrogram(hclust_ward)

# ggdend
gdend_ward <- as.ggdend(dend_ward)
gdend_ward$labels$angle <- 90

# add colors
gdend_ward$labels %>% left_join(sarc_GSM_to_col_final, by ="label") %>%
  select(1,2,5,6,7,8) %>%
  rename("label" = 5, "col" = 6) ->gdend_ward$labels

save(sarc_GSM_to_col_final, hclust_ward, dend_ward, gdend_ward, file = file.path("results", "sarc_ref_val_HC_ward.RData"))

# draw png
pngpath <- "results/sarc_hclust_ward.png"
png(file = pngpath, width = 19000, height = 10000)
ggplot(gdend_ward)
dev.off()

############################################
# method = complete and circlize
gdend_ward_circle <- gdend_ward
gdend_ward_circle$labels <- gdend_ward$labels %>% mutate(angle = (90 - 360 * x /(max(x) + 1)))


pngpath <- "results/sarc_hclust_ward_circ.png"
png(file = pngpath, width = 20000, height = 20000)
ggplot() +
  geom_segment(data = gdend_ward$segments, aes(x =x, y = y, xend = xend, yend = yend, col = col)) +
  geom_text(data = gdend_ward_circle$labels, aes(x =x, y = y-10, label = label, col = col, angle = angle),hjust = 0) +
  scale_x_continuous(expand = c(0,1))+
  scale_y_reverse(expand = c(0,100)) +
  coord_polar(theta="x") +
  theme_minimal()+ 
  theme(legend.position = "none",
        axis.line.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.title.y=element_blank()) + 
  scale_color_identity()
dev.off()

####################################################
# add black

##########################################
# ggdend
gdend_comp_bl <- as.ggdend(dend_comp)
gdend_comp_bl$labels$angle <- 90


# make table GSM to color black!
sarc_GSM_to_col_final_bl <- sarc_GSM_to_col_final
sarc_GSM_to_col_final_bl[1:1077,2] ="#000000"

# add colors
gdend_comp_bl$labels %>% left_join(sarc_GSM_to_col_final_bl, by ="label") %>%
  select(1,2,5,6,7,8) %>%
  rename("label" = 5, "col" = 6) ->gdend_comp_bl$labels


# draw png
pngpath <- "results/sarc_hclust_comp_bl.png"
png(file = pngpath, width = 19000, height = 10000)
ggplot(gdend_comp_bl)
dev.off()


#########################################
# method = complete and circlize
gdend_comp_bl_circle <- gdend_comp_bl
gdend_comp_bl_circle$labels <- gdend_comp_bl$labels %>% mutate(angle = (90 - 360 * x /(max(x) + 1)))


pngpath <- "results/sarc_hclust_comp_bl_circ.png"
png(file = pngpath, width = 20000, height = 20000)
ggplot() +
  geom_segment(data = gdend_comp_bl$segments, aes(x =x, y = y, xend = xend, yend = yend, col = col)) +
  geom_text(data = gdend_comp_bl_circle$labels, aes(x =x, y = y-10, label = label, col = col, angle = angle),hjust = 0) +
  scale_x_continuous(expand = c(0,1))+
  scale_y_reverse(expand = c(0,100)) +
  coord_polar(theta="x") +
  theme_minimal()+ 
  theme(legend.position = "none",
        axis.line.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.title.y=element_blank()) + 
  scale_color_identity()
dev.off()


#########################################


###############################################
# method = ward.d2

# ggdend
gdend_ward_bl <- as.ggdend(dend_ward)
gdend_ward_bl$labels$angle <- 90

# add colors
gdend_ward_bl$labels %>% left_join(sarc_GSM_to_col_final_bl, by ="label") %>%
  select(1,2,5,6,7,8) %>%
  rename("label" = 5, "col" = 6) ->gdend_ward_bl$labels

# draw png
pngpath <- "results/sarc_hclust_ward_bl.png"
png(file = pngpath, width = 19000, height = 10000)
ggplot(gdend_ward_bl)
dev.off()

############################################
# method = complete and circlize
gdend_ward_bl_circle <- gdend_ward_bl
gdend_ward_bl_circle$labels <- gdend_ward_bl$labels %>% mutate(angle = (90 - 360 * x /(max(x) + 1)))

pngpath <- "results/sarc_hclust_ward_bl_circ.png"
png(file = pngpath, width = 20000, height = 20000)
ggplot() +
  geom_segment(data = gdend_ward_bl$segments, aes(x =x, y = y, xend = xend, yend = yend, col = col)) +
  geom_text(data = gdend_ward_bl_circle$labels, aes(x =x, y = y-10, label = label, col = col, angle = angle),hjust = 0) +
  scale_x_continuous(expand = c(0,1))+
  scale_y_reverse(expand = c(0,100)) +
  coord_polar(theta="x") +
  theme_minimal()+ 
  theme(legend.position = "none",
        axis.line.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank(),
        axis.title.y=element_blank()) + 
  scale_color_identity()
dev.off()

####################################################


