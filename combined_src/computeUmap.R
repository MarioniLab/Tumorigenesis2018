# Compute UMAP for whole dataste
library(dplyr)
library(ggplot2)
library(Matrix)
library(umap)
library(scran)
source("functions.R")

dataList <- readRDS("../data/combined_Robjects/ExpressionList_QC_norm.rds")
# dataList <- subSample(dataList,cell.number=60000)
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
rm(dataList)

# Highly variable genes
hvg <- getHighVar(m,get.var.out=TRUE, supress.plot=TRUE)
hvg <- hvg[rownames(hvg) %in% fD$uniqnames[fD$KeepForHvg],]
hvg <- rownames(hvg[order(hvg$bio, decreasing=TRUE),])[1:(nrow(hvg)/10)]

ump <- umap(as.matrix(t(m[hvg,])), min_dist=0.5, metric="pearson", random_state=42)
pD$UMAP1Uncor <- ump$layout[,1]
pD$UMAP2Uncor <- ump$layout[,2]

m.cor <- readRDS("../data/combined_Robjects/CorrectedPCA.rds")
m.cor <- m.cor[pD$barcode,]
ump <- umap(m.cor, min_dist=0.5, metric="pearson", random_state=42)
pD$UMAP1 <- ump$layout[,1]
pD$UMAP2 <- ump$layout[,2]

write.csv(pD[,c("barcode","UMAP1","UMAP2","UMAP1Uncor","UMAP2Uncor")],file="../data/combined_Robjects/UMAP_corrected.csv")
