# Compute UMAP for whole dataste
library(dplyr)
library(ggplot2)
library(Matrix)
library(umap)
library(scran)
source("functions.R")

dataList <- readRDS("../data/firstBatch_Robjects/ExpressionList_QC.rds")
sfs <- read.csv("../data/firstBatch_Robjects/SizeFactors.csv")
# dataList <- subSample(dataList,cell.number=10000)
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
pD <- right_join(pD,sfs)
rm(dataList)

m <- m[,pD$barcode]
m <- log2(t(t(m)/pD$sf)+1)
m <- as(m,"dgCMatrix")

fD <- fD[Matrix::rowMeans(m)>0.01,] 
m <- m[Matrix::rowMeans(m)>0.01,]

# Highly variable genes
hvg <- getHighVar(m,get.var.out=TRUE, supress.plot=TRUE)
hvg <- hvg[rownames(hvg) %in% fD$uniqnames[fD$KeepForHvg],]
hvg <- rownames(hvg[order(hvg$bio, decreasing=TRUE),])[1:(nrow(hvg)/10)]

ump <- umap(as.matrix(t(m[hvg,])), min_dist=0.5, metric="pearson", random_state=42)
pD$UMAP1 <- ump$layout[,1]
pD$UMAP2 <- ump$layout[,2]

write.csv(pD[,c("barcode","UMAP1","UMAP2")],file="../data/firstBatch_Robjects/UMAP_all.csv")
