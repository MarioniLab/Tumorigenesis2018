# Load Data
library(dplyr)
library(igraph)
library(scran)
library(umap)
library(Matrix)
source("functions.R")
dataList <- readRDS("../data/combined_Robjects/ExpressionList_QC_norm.rds")
pD <- dataList[["phenoData"]]
m <- dataList[["counts"]]
fD <- dataList[["featureData"]]
rownames(pD) <- pD$barcode
rm(dataList)
ump <- read.csv("../data/combined_Robjects/UMAP_corrected.csv",stringsAsFactors=FALSE)
pD <- pD[ump$barcode,] # this ensure the right order and cells for the umap graph below
m <- m[,pD$barcode]

# ---- First Round -----
# Clustering with walktrap on the UMAP graph
umapgraph <- readRDS("../data/combined_Robjects/UMAP_graphs.rds")
igr.cor <- umapgraph[["Corrected"]]
set.seed(42)
cl.cor <- cluster_walktrap(igr.cor,steps=5)
pD$Cluster <- as.factor(paste0("C",cl.cor$membership))
write.csv(file="../data/combined_Robjects/Clusters.csv",pD[,c("barcode","Cluster")])

## Subclustering clusters

mergeCluster <- function(x, clusters, min.DE=10, maxRep=30, removeGenes=NULL, merge=TRUE, ...)
{
    library(scran)
    counter <- c(1:maxRep)
    if (!is.null(removeGenes)) {
	x <- x[!(rownames(x) %in% removeGenes),]
    }
    for (i in counter) {
	clust.vals <- levels(as.factor(clusters))
	marker.list <- findMarkers(x,clusters,subset.row=rowMeans(x)>0.01, lfc=1, full.stats=TRUE, ...)

	out <- matrix(nrow=length(clust.vals), ncol=length(clust.vals))
	colnames(out) <- clust.vals
	rownames(out) <- clust.vals
	for (cl in names(marker.list)) {
	    sb <- marker.list[[cl]]
	    sb <- data.frame(sb)
	    sb <- as.data.frame(sb[,grepl(".log.FDR",colnames(sb))]) # as.data.frame in case it's a single column
	    if(ncol(sb)==1) {
		colnames(sb) <- setdiff(names(marker.list),cl)
	    } else {
	    colnames(sb) <- gsub("stats.|.log.FDR","",colnames(sb))
	    }
	    sb <- colSums(as.matrix(sb) < -2)
	    sb[cl] <- 0
	    sb <- sb[match(colnames(out),names(sb))]
	    out[,cl] <- sb
	}

	if (merge) {
	#Compute new groups
	if (any(is.na(out))) { print("A test did not have enough D.F.") }
	out[is.na(out)] <- min.DE+1 # This is to prevent hclust from crashing, NAs are produced when their are not enough d.f. for testing. In this case I do not wnat the clusters to be merged because I formally cant test for DE.
	hc <- hclust(as.dist(out))
	min.height <- min(hc$height)
	if (min.height <= min.DE) {
	    print(paste0(i,". joining of clusters with ",min.height," DE genes"))
	    cs <- cutree(hc,h=min.height)
	    newgrps <- sapply(c(1:max(cs)), function(i) paste(names(cs)[cs==i],collapse="."))
	    csnew <- plyr::mapvalues(cs,c(1:max(cs)), newgrps)
	    clusters <- plyr::mapvalues(clusters,names(csnew),csnew)
	} else {
	    break
	}
    } else {
	break}

	if (length(unique(clusters))==1) {
	    break
	}
    }
	output <- list("Marker"=marker.list,
		       "NumberDE"=out,
		       "NewCluster"=clusters)
	return(output)
}


compCluster <-  function(pcs, cluster, m, ...) {

	if (nrow(pcs)<100) {
	    out <- data.frame("barcode"=rownames(pcs),
			      "MajorCluster"=cluster,
			      "Cluster"=cluster,
			      "SCID"=1,
			      "UnmergedID"=1,
			      "SubUMAP1"=NA,
			      "SubUMAP2"=NA)
	} else {

	    ump <- umap(pcs, random_state=42)
	    umap1 <- ump$layout[,1]
	    umap2 <- ump$layout[,2]

	    igr <- get_umap_graph(ump)

	    set.seed(42)
	    wktrp <- cluster_walktrap(igr)
	    subclusters <- wktrp$membership

	    # Merging step
	    if (length(unique(subclusters))>1) {
		mrg <- mergeCluster(m, subclusters, ...)
		finalSubclusters <- as.character(mrg$NewCluster)
		finalSubclusters <- plyr::mapvalues(finalSubclusters,unique(finalSubclusters),c(1:length(unique(finalSubclusters))))
	    } else {
		finalSubclusters <- subclusters
	    }

	    cs <- paste0(cluster,"S", finalSubclusters)

	    out <- data.frame("barcode"=rownames(pcs),
			      "MajorCluster"=cluster,
			      "Cluster"=cs,
			      "SCID"=finalSubclusters,
			      "UnmergedID"=subclusters,
			      "SubUMAP1"=umap1,
			      "SubUMAP2"=umap2)
	}
	return(out)
}



# Load Data
cor.pca <- readRDS("../data/combined_Robjects/CorrectedPCA.rds")
cor.pca <- cor.pca[pD$barcode,]
rownames(pD) <- pD$barcode

clusts <- as.character(unique(pD$Cluster))
## remove later
out <- bplapply(clusts, function(cl) {
	     print(cl)
	     cells <- pD$barcode[pD$Cluster==cl]
	     pD.sub <- pD[cells,]
	     cor.pca.sub <- cor.pca[cells,]
	     m.sub <- m[,cells]
	     out <- compCluster(cor.pca.sub, cluster=cl, m=m.sub, block=pD.sub$Batch, removeGenes=rownames(fD)[!fD$KeepForHvg])
	     return(out)
			  }, BPPARAM=MulticoreParam(4))

out.tbl <- do.call(rbind,out)
out.tbl$Cluster <- as.character(out.tbl$Cluster)
out.tbl$Cluster <- factor(out.tbl$Cluster,levels=gtools::mixedsort(unique(out.tbl$Cluster)))
write.csv(file="../data/combined_Robjects/SubClusters.csv",out.tbl)
