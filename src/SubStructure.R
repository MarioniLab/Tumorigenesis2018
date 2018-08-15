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


dataList <- readRDS("../data/Robjects/ExpressionList_QC.rds")
# dataList <- subSample(dataList,cell.number=10000)
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
rm(dataList)
sfs <- read.csv("../data/Robjects/SizeFactors.csv")
cluster <- read.csv("../data/Robjects/CellTypes.csv")

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

compTsne <-  function(m, fD) {
	# Highly variable genes
	hvg <- getHighVar(m,get.var.out=TRUE, supress.plot=TRUE)
	hvg <- hvg[rownames(hvg) %in% fD$id[fD$KeepForHvg],]
	hvg <- rownames(hvg[order(hvg$bio, decreasing=TRUE),])[1:(nrow(hvg)/10)]

	# Compute tSNE 
	fPCA <- m[hvg,]
	set.seed(300)
	pca <- prcomp_irlba(t(fPCA), n=50, center=FALSE, scale.=FALSE)
	pca <- pca$x
	set.seed(300)
	tsn <- Rtsne(pca,perplexity=50,pca=FALSE)
	out <- data.frame("barcode"=colnames(m),
			  "tSNE1"=tsn$Y[,1],
			  "tSNE2"=tsn$Y[,2])
	return(out)
}


conds <- list("Epithelial","Fibroblast",c("Lymphoid","Myeloid"))
out <- data.frame()

for (cond in conds) {
    m.sub <- m[,pD$barcode[pD$Class %in% cond]]
    tmp <- compTsne(m.sub, fD)
    out <- rbind(out,tmp)
}

pD.sub <- right_join(pD,out)

pD.sub$Class <- as.character(pD.sub$Class)
pD.sub$Class[pD.sub$Class %in% c("Lymphoid","Myeloid")] <- "Immune"
pD.sub$Class <- factor(pD.sub$Class)

res <- pD.sub[,c("barcode","Cluster","Class","tSNE1","tSNE2","Age")]
write.csv(file="../data/Robjects/Substructure.csv",res)
