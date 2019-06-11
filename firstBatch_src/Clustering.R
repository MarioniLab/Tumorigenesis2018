# First level Clustering using SNN Graph
library(RColorBrewer)
library(reshape2)
library(plyr)
library(dplyr)
library(scran)
library(Rtsne)
library(ggplot2)
library(Matrix)
library(igraph)
library(scater)
source("functions.R")


# Load Data
dataList <- readRDS("../data/firstBatch_Robjects/ExpressionList_QC.rds")
# dataList <- subSample(dataList,cell.number=10000)
m.full <- dataList[["counts"]]
pD.full <- dataList[["phenoData"]]
fD.full <- dataList[["featureData"]]
rm(dataList)
sfs <- read.csv("../data/firstBatch_Robjects/SizeFactors.csv")
rownames(m.full) <- as.vector(rownames(m.full))

# Gene and cell filtering
m <- m.full[,pD.full$PassAll]
keep <- rowMeans(m)>0.01
pD <- pD.full[pD.full$PassAll,]
m <- m[keep,]
fD <- fD.full[keep,]

#Normalize
pD <- left_join(pD,sfs)
m <- log2(t(t(m)/pD$sf)+1)
m <- as(m,"dgCMatrix")

# Highly Variable Genes
hvg <- getHighVar(m,get.var.out=TRUE,supress.plot=TRUE)
hvg <- hvg[(rownames(hvg) %in% fD$uniqnames[fD$KeepForHvg]),]	
hvg <- rownames(hvg[order(hvg$bio, decreasing=TRUE),])[1:(nrow(hvg)/10)]

# Clustering
set.seed(300)
igr <- buildSNNGraph(m[hvg,],pc.approx=TRUE,k=10)
cs <- cluster_louvain(igr)
cs <- cs$membership
clusts <- as.factor(paste0("C",unname(cs)))
pD$Cluster <- clusts

# Save
anno <- pD[,c("barcode","Cluster")]
write.csv(file="../data/firstBatch_Robjects/Cluster_all.csv",anno)

