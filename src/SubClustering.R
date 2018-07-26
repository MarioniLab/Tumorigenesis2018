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
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
rm(dataList)
sfs <- read.csv("../data/Robjects/SizeFactors.csv")
cluster <- read.csv("../data/Robjects/Cluster_all.csv")

rownames(m) <- as.vector(rownames(m))

pD <- left_join(pD,cluster)
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

	# Compute Clusters 

	igr <- buildSNNGraph(pca,pc.approx=TRUE,d=NA,transpose=TRUE)
	cs <- cluster_louvain(igr)
	cs <- cs$membership
	clusts <- as.factor(paste0("SC",unname(cs)))
	merged <- mergeCluster(m, clusts, maxRep=40, min.DE=15, removeGenes=fD$symbol[!fD$KeepForHvg])
	cluster <- merged$NewCluster

	out <- data.frame("barcode"=colnames(m),
			  "Sub-tSNE1"=tsn$Y[,1],
			  "Sub-tSNE2"=tsn$Y[,2],
			  "Cluster"=cluster)
	return(out)
}

clsters <- levels(pD$Cluster)
out <- data.frame()

for (clsters in clusters) {
    m.sub <- m[,pD$barcode[pD$Cluster %in% cluster]]
    tmp <- compTsne(m.sub, fD)
    out <- rbind(out,tmp)
    print(paste0("Done with", cluster))
}

pD.sub <- right_join(pD,out)
