library(RColorBrewer)
library(reshape2)
library(plyr)
library(dplyr)
library(scran)
library(Rtsne)
library(ggplot2)
library(Matrix)
library(igraph)
library(dynamicTreeCut)
library(cluster)
library(cowplot)
library(scater)
library(irlba)
source("functions.R")


dataList <- readRDS("../data/firstBatch_Robjects/ExpressionList_QC.rds")
# dataList <- subSample(dataList,cell.number=10000)
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
rm(dataList)
sfs <- read.csv("../data/firstBatch_Robjects/SizeFactors.csv")
cluster <- read.csv("../data/firstBatch_Robjects/Cluster_all.csv")

rownames(m) <- as.vector(rownames(m))

pD <- right_join(pD,cluster)
pD <- left_join(pD,sfs)

# Gene and cell filtering
m <- m[,pD$barcode]
keep <- rowMeans(m)>0.01
m <- m[keep,]
fD <- fD[keep,]

# Norm and Log
m <- log2(t(t(m)/pD$sf)+1)

cells <- pD$barcode[pD$Class=="Epithelial" & !(pD$Condition %in% c("CBLT","WKBR61.4dT"))]
epi <- pD[pD$barcode %in% cells,]
m <- m[,epi$barcode]

fD <- fD[rowMeans(m)>0.01,]
m <- m[rowMeans(m)>0.01,]

# 
hvg <- getHighVar(m,get.var.out=TRUE,supress.plot=TRUE)
hvg <- hvg[(rownames(hvg) %in% fD$id[fD$KeepForHvg]),]	
hvg <- rownames(hvg[order(hvg$bio, decreasing=TRUE),])[1:(nrow(hvg)/10)]

#Diffmap
library(destiny)
fDmap <- t(as.matrix(m[hvg,]))
set.seed(300)
dm <- DiffusionMap(fDmap)#,rotate=TRUE)
#test 
cluster  <-  destiny:::get_louvain_clusters(dm@transitions)
epi$DMcluster <- factor(cluster)
dpt <- DPT(dm)

dms <- eigenvectors(dm)[,1:5]
dms <- data.frame(dms,
		  barcode=epi$barcode,
		  branch=dpt$branch)
epi <- left_join(epi,dms)
ggplot(epi, aes(x=DC1,y=DC2,color=branch)) +
    geom_point() 

out <- epi[,c("barcode","DC1","DC2","DC3","DC4","DC5","branch","DMcluster")]
write.csv(out,file="../data/firstBatch_Robjects/DiffusionMap_Epithelium_Cluster.csv")

# subs <- epi$barcode[epi$DMcluster %in% c(2,3,5)]
# m.sub <- m[,subs]
# epi.sub <- epi[epi$barcode %in% subs,]
# markers <- findMarkers(m.sub,factor(epi.sub$DMcluster),direction="up",lfc=log2(1.3))
# x <-  markers[[3]]
# x$id <- rownames(x)
# x <- left_join(data.frame(x),fD[,c("id","symbol")])

# sce <- SingleCellExperiment(rowData=fD,
#                             colData=epi,
#                             assays=list(logcounts=m),
#                             reducedDims=SimpleList(DC1=as.matrix(epi[,c("DC1","DC2")]),
#                                                    DC2=as.matrix(epi[,c("DC2","DC3")])))
# require(iSEE)
# iSEE(sce)
# 
