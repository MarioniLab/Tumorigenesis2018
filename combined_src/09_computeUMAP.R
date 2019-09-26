# Compute UMAP for whole dataste
library(dplyr)
library(ggplot2)
library(Matrix)
library(umap)
library(scran)
library(igraph)
source("functions.R")

dataList <- readRDS("../data/combined_Robjects/ExpressionList_QC_norm.rds")
# dataList <- subSample(dataList,cell.number=60000)
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
rm(dataList)

# Remove Doublets 
qcMets <- read.csv("../data/combined_Robjects/QC_Part2.csv",row.names=1,stringsAsFactors=FALSE)
rmCells <- qcMets$barcode[qcMets$isDoubletFinal | qcMets$isRbc]

pD <- pD[!pD$barcode %in% rmCells,]
m <- m[,pD$barcode]

# ---- Pre-Correction ----

# Highly variable genes
hvg <- getHighVar(m,get.var.out=TRUE, supress.plot=TRUE)
hvg <- hvg[rownames(hvg) %in% fD$uniqnames[fD$KeepForHvg],]
hvg <- rownames(hvg[order(hvg$bio, decreasing=TRUE),])[1:(nrow(hvg)/10)]

# UMAP
ump <- umap(as.matrix(t(m[hvg,])), random_state=42)
pD$UMAP1Uncor <- ump$layout[,1]
pD$UMAP2Uncor <- ump$layout[,2]

#Convert to adj matrix for clustering if desired
indx.knn <- ump$knn[["indexes"]]
m.adj <- Matrix(0, nrow(indx.knn), nrow(indx.knn)) 
rownames(m.adj) <- colnames(m.adj) <- rownames(indx.knn)

for (i in seq_len(nrow(m.adj))) {
    m.adj[i,rownames(indx.knn)[indx.knn[i,]]] <- 1
}

igr.uncor <- graph_from_adjacency_matrix(m.adj, mode="undirected")#,weighted=TRUE)

# ---- After Correction -----
m.cor <- readRDS("../data/combined_Robjects/CorrectedPCA.rds")
m.cor <- m.cor[pD$barcode,]
ump <- umap(m.cor, random_state=42)
pD$UMAP1 <- ump$layout[,1]
pD$UMAP2 <- ump$layout[,2]

#Convert to adj matrix for clustering if desired
indx.knn <- ump$knn[["indexes"]]
m.adj <- Matrix(0, nrow(indx.knn), nrow(indx.knn)) 
rownames(m.adj) <- colnames(m.adj) <- rownames(indx.knn)

for (i in seq_len(nrow(m.adj))) {
    m.adj[i,rownames(indx.knn)[indx.knn[i,]]] <- 1
}

igr.cor <- graph_from_adjacency_matrix(m.adj, mode="undirected")#,weighted=TRUE)


out <- list("Corrected"=igr.cor,
	    "Uncorrected"=igr.uncor)
saveRDS(out,"../data/combined_Robjects/UMAP_graphs.rds")
write.csv(pD[,c("barcode","UMAP1","UMAP2","UMAP1Uncor","UMAP2Uncor")],file="../data/combined_Robjects/UMAP_corrected.csv")
