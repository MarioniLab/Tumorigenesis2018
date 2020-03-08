# Clustering
library(igraph)
library(scran)
library(Matrix)
source("functions.R")

# Load Data
sce <- readRDS("../data/combined_Robjects/SCE_QC_norm.rds")
ump <- read.csv("../data/combined_Robjects/UMAP_corrected.csv",stringsAsFactors=FALSE)[,-1]
sce <- sce[,ump$barcode] # this ensure the right order and cells for the umap graph below

# Graph
umapgraph <- readRDS("../data/combined_Robjects/UMAP_graphs.rds")
igr.cor <- umapgraph[["Corrected"]]

# ---- Based on UMAP Graph ----
# Clustering with walktrap on the UMAP graph
set.seed(42)
ump.wktrp <- cluster_walktrap(igr.cor, steps=5)
cl.ump.wktrp <- paste0("C",ump.wktrp$membership)


# ---- Merge Clusters ----
rmgenes <- rownames(sce)[!rowData(sce)$KeepForHvg]
sce$Cluster <- cl.ump.wktrp
m <- logcounts(sce)
merged <- mergeCluster(m, factor(sce$Cluster), removeGenes=rmgenes, block=sce$Batch)
sce$Cluster <- merged$NewCluster

# ---- Save ----
out <- colData(sce)
out$Cluster <- paste0("C",as.numeric(out$Cluster))
write.csv(out[,c("barcode","Cluster")],"../data/combined_Robjects/Clusters.csv")
