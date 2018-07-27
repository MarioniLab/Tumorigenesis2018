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

pD <- right_join(pD,cluster) # This also "performs" QC
pD <- left_join(pD,sfs)

# Gene and cell filtering
m <- m[,pD$barcode]
keep <- rowMeans(m)>0.01
m <- m[keep,]
fD <- fD[keep,]

# Norm and Log
m <- log2(t(t(m)/pD$sf)+1)

compTsne <-  function(m, fD, cluster,clustering="graph") {
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
	tsn <- Rtsne(pca,perplexity=30,pca=FALSE)

	# Compute Clusters 

	if (clustering=="graph") {
	    igr <- buildSNNGraph(pca,pc.approx=TRUE,d=NA,transpose=TRUE, k=30)
	    cs <- cluster_louvain(igr)
	    cs <- cs$membership
	    clusts <- as.factor(paste0("SC",unname(cs)))
	    merged <- mergeCluster(m, clusts, maxRep=40, min.DE=15, removeGenes=fD$symbol[!fD$KeepForHvg])
	} else {
	    dis <- dist(pca)
	    tree <- hclust(dis, method="ward.D2")
	    maxCoreScatter <- 0.3
	    minGap <- (1 - maxCoreScatter) * 3/4
	    cs <- cutreeDynamic(tree,distM=as.matrix(dis),deepSplit=0, minClusterSize=max(50,floor(nrow(pca)/100)), maxCoreScatter=maxCoreScatter, minGap=minGap)
	    cs  <- as.factor(paste0("SC",unname(cs)))
	    merged <- mergeCluster(m, cs, maxRep=40, min.DE=15, removeGenes=fD$symbol[!fD$KeepForHvg])
	}

	cs <- paste0(cluster, merged$NewCluster)
	out <- data.frame("barcode"=colnames(m),
			  "Sub-tSNE1"=tsn$Y[,1],
			  "Sub-tSNE2"=tsn$Y[,2],
			  "Cluster"=cs)
	return(out)
}

clsters <- levels(pD$Cluster)
clsters <- clsters[clsters!="C8"]
out <- data.frame()

for (cluster in clsters) {
    m.sub <- m[,pD$barcode[pD$Cluster %in% cluster]]
    tmp <- compTsne(m.sub, fD, cluster)#, clustering ="hierarchical")
    out <- rbind(out,tmp)
    print(paste0("Done with", cluster))
}

colnames(out)[4] <- "SubCluster"
pD.sub <- right_join(pD,out,by="barcode")

plotlist <- list()
for (clster in levels(pD.sub$Cluster)) {
    pD.sub.sub <- pD.sub[pD.sub$Cluster==clster,]
    plotlist[[clster]] <- ggplot(pD.sub.sub, aes(x=Sub.tSNE1, y=Sub.tSNE2, color=SubCluster)) +
	geom_point () +
	theme_void()
}

pD.sub$SubClusterLabels <- mapvalues(pD.sub$SubCluster,
				 levels(pD.sub$SubCluster),
                         	 c("C1A","C1B","C1C","C1D","C1E",
				   "C10A","C10B","C11A","C11B","C11C",
				   "C12","C13A","C13B","C13C","C13D",
				   "C14A","C15","C16A","C16B","C16C","C16D",
				   "C17","C18A","C18B","C18C","C18D","C18E",
				   "C18F","C18G","C2A","C2B","C2C","C19A","C19B",
				   "C19C","C19D","C20A","C20B","C3A","C3B",
				   "C3C","C3D","C4A","C4B","C4C","C4D","C4E",
				   "C5A","C5B","C5C","C5D","C6A","C6B","C6C",
				   "C7","C9A","C9B","C9C","C9D"))
out <- pD.sub[,c("barcode","SubClusterLabels","Sub.tSNE1","Sub.tSNE2")]
res <- left_join(pD,out)
res$SubClusterLabels <- as.character(res$SubClusterLabels)
res$SubClusterLabels[is.na(res$SubClusterLabels)] <- "C8"
res$ClusterLabels <- gsub("A$|B$|C$|D$|E$|F$|G$|H$","",res$SubClusterLabels)

res <- res[,c("barcode","ClusterLabels","SubClusterLabels","Sub.tSNE1","Sub.tSNE2")]
write.csv(res,"../data/Robjects/Cluster_final.csv")
