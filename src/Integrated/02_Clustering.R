# Clustering
library(igraph)
library(scran)
library(Matrix)
source("../functions.R")

# Load Data
sce <- readRDS("../../data/Integrated/Robjects/Pregnancy_FullMNN.rds")

# Graph
igr.cor <- readRDS("../../data/Integrated/Robjects/PregnancyUMAP_graph.rds")

# ---- Based on UMAP Graph ----
# Clustering with walktrap on the UMAP graph
set.seed(42)
ump.wktrp <- cluster_walktrap(igr.cor, steps=7)
cl.ump.wktrp <- paste0("C",ump.wktrp$membership)


# ---- Merge Clusters ----
rD <- rowData(readRDS("../../data/Tumorigenesis/Robjects/SCE_QC_norm.rds"))
rmgenes <- rownames(rD)[!rD$KeepForHvg]
rmgenes <- rownames(sce) %in% rmgenes
sce$Cluster <- cl.ump.wktrp
m <- logcounts(sce)
merged <- mergeCluster(m, factor(sce$Cluster), removeGenes=rmgenes, block=sce$Batch)
sce$Cluster <- merged$NewCluster

# ---- Save ----
out <- colData(sce)
out$Cluster <- paste0("C",as.numeric(out$Cluster))
write.csv(out[,c("barcode","Cluster")],"../../data/Integrated/Robjects/PregnancyClusters.csv")
