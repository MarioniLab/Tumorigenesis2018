# SubClustering
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


# Load Data
dataList <- readRDS("../data/firstBatch_Robjects/ExpressionList_QC.rds")
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
rm(dataList)
sfs <- read.csv("../data/firstBatch_Robjects/SizeFactors.csv")
cluster <- read.csv("../data/firstBatch_Robjects/Cluster_all.csv")

rownames(m) <- as.vector(rownames(m))

pD <- right_join(pD,cluster) # This also "performs" QC
pD <- left_join(pD,sfs)

# Gene and cell filtering
m <- m[,pD$barcode]
keep <- rowMeans(m)>0.01
m <- m[keep,]
fD <- fD[keep,]

# Norm and Log
m <- log2(t(t(m)/pD$sf)+1)
m <- as(m, "dgCMatrix")

compCluster <-  function(m, fD, cluster,clustering="graph") {
	# Highly variable genes
	hvg <- getHighVar(m,get.var.out=TRUE, supress.plot=TRUE)
	hvg <- hvg[rownames(hvg) %in% fD$uniqnames[fD$KeepForHvg],]
	hvg <- rownames(hvg[order(hvg$bio, decreasing=TRUE),])[1:(nrow(hvg)/10)]

	# Compute PCA
	fPCA <- m[hvg,]
	set.seed(300)
	pca <- prcomp_irlba(t(fPCA), n=50, center=FALSE, scale.=FALSE)
	pca <- pca$x

	# Compute Clusters 

	if (clustering=="graph") {
	    igr <- buildSNNGraph(pca, d=NA,transpose=TRUE, k=20)
	    cs <- cluster_louvain(igr)
	    cs <- cs$membership
	    clusts <- as.factor(paste0("SC",unname(cs)))
	    #             merged <- mergeCluster(m, clusts, maxRep=40, min.DE=10, removeGenes=fD$uniqnames[!fD$KeepForHvg])
	} else {
	    dis <- dist(pca)
	    tree <- hclust(dis, method="ward.D2")
	    maxCoreScatter <- 0.3
	    minGap <- (1 - maxCoreScatter) * 3/4
	    cs <- cutreeDynamic(tree,distM=as.matrix(dis),deepSplit=0, minClusterSize=max(50,floor(nrow(pca)/100)), maxCoreScatter=maxCoreScatter, minGap=minGap)
	    cs  <- as.factor(paste0("SC",unname(cs)))
	    merged <- mergeCluster(m, cs, maxRep=40, min.DE=15, removeGenes=fD$symbol[!fD$KeepForHvg])
	}

	#         cs <- paste0(cluster, merged$NewCluster)
	cs <- paste0(cluster, clusts)
	out <- data.frame("barcode"=colnames(m),
			  "Cluster"=cs)
	return(out)
}

# Exclude C8 because too small
clsters <- levels(pD$Cluster)
clsters <- clsters[clsters!="C8"]
out <- data.frame()

for (cluster in clsters) {
    m.sub <- m[,pD$barcode[pD$Cluster %in% cluster]]
    tmp <- compCluster(m.sub, fD, cluster)
    out <- rbind(out,tmp)
    print(paste0("Done with", cluster))
}

# Reformat
colnames(out)[4] <- "SubCluster"
pD.sub <- right_join(pD,out,by="barcode")
pD <- left_join(pD,out,by="barcode")
pD$SubCluster <- as.character(pD$SubCluster)
pD$SubCluster[is.na(pD$SubCluster)] <- "C8"
pD$SubCluster <- as.factor(pD$SubCluster)

# Merge Clusters with less than 15 DE Genes
merged <- mergeCluster(m, pD$SubCluster, maxRep=40, min.DE=15, removeGenes=fD$uniqnames[!fD$KeepForHvg])
pD$NewCluster <- merged$NewCluster

# Rename in a sensible way
library(gtools)
pD$FinalCluster <- paste0("C",as.numeric(factor(pD$NewCluster,levels=mixedsort(levels(pD$NewCluster)))))

# Save
res <- pD[,c("barcode","Cluster","SubCluster","NewCluster","FinalCluster")]
write.csv(res,"../data/firstBatch_Robjects/Cluster_final.csv")
