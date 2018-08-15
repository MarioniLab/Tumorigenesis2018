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
source("functions.R")


dataList <- readRDS("../data/Robjects/ExpressionList_QC.rds")
# dataList <- subSample(dataList,cell.number=10000)
m.full <- dataList[["counts"]]
pD.full <- dataList[["phenoData"]]
fD.full <- dataList[["featureData"]]
rm(dataList)
sfs <- read.csv("../data/Robjects/SizeFactors.csv")

rownames(m.full) <- as.vector(rownames(m.full))
# cond <- c("38","41","53")

cells <- pD.full$Condition!="WKBR61.4dT"
# Gene and cell filtering
m <- m.full[,pD.full$PassAll & cells]
keep <- rowMeans(m)>0.01
pD <- pD.full[pD.full$PassAll & cells,]
m <- m[keep,]
fD <- fD.full[keep,]

#Normalize
pD <- left_join(pD,sfs)

m <- log2(t(t(m)/pD$sf)+1)

hvg <- getHighVar(m,get.var.out=TRUE,supress.plot=TRUE)
hvg <- hvg[(rownames(hvg) %in% fD$id[fD$KeepForHvg]),]	
hvg <- rownames(hvg[order(hvg$bio, decreasing=TRUE),])[1:(nrow(hvg)/10)]

# clustering
set.seed(300)
igr <- buildSNNGraph(m[hvg,],pc.approx=TRUE)
cs <- cluster_louvain(igr)
cs <- cs$membership
clusts <- as.factor(paste0("C",unname(cs)))

# Merge Clusters with less DE Genes than 15
merged <- mergeCluster(m, clusts, maxRep=40, min.DE=15, removeGenes=fD$id[!fD$KeepForHvg])
pD$Cluster <- merged$NewCluster

# Test Graph abstraction
x <- clusterModularity(igr,pD$Cluster,get.values=TRUE)
obexp <- x$observed/x$expected
obexp[obexp <= 10^-1] <- 0
grph <- graph_from_adjacency_matrix(obexp,weighted=TRUE,diag=FALSE,mode="undirected")
plot(grph,edge.width=E(grph)$weight*20, vertex.size=as.vector(table(pD$Cluster))/300)

# Compute tSNE 
tsn <- read.csv("../data/Robjects/tSNE_nontumor_raw.csv")
pD <- left_join(pD,tsn)


p1 <- ggplot(pD, aes(x=tSNE1, y=tSNE2, color=Cluster)) +
    geom_point() +
    theme_void() +
    facet_wrap(~Cluster)
p1

anno <- pD[,c("barcode","Cluster")]
anno$Class <- mapvalues(anno$Cluster,
			c("C1","C2.C3","C4","C5","C6","C7","C8","C9","C10",
			  "C11","C12","C13","C14","C15","C16","C17","C18"),
			c("Lymphoid","Lymphoid","Lymphoid","Fibroblast","Lymphoid","Fibroblast","Epithelial","Fibroblast","Epithelial",
			  "Myeloid","Lymphoid","Lymphoid","Myeloid","Fibroblast","Epithelial","Doublet","Endothelial"))


write.csv(file="../data/Robjects/CellTypes.csv",anno)

# rownames(pD) <- pD$barcode
# rownames(fD) <- fD$symbol
# rownames(m) <- fD$symbol
# marker <- findMarkers(m,pD$Cluster)
# sapply(names(marker),function(x) write.csv(marker[[x]],file=paste0("tmp/",x,".csv")))

# sce <- SingleCellExperiment(rowData=fD,
#                             colData=pD,
#                             assays=list(logcounts=m),
#                             reducedDims=SimpleList(tSNE=as.matrix(pD[,c("tSNE1","tSNE2")])))
# 
# library(iSEE)
# 
# iSEE(sce)
# 
