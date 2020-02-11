# Clustering
# library(dplyr)
library(igraph)
library(scran)
# library(umap)
library(Matrix)
source("functions.R")

# Load Data
sce <- readRDS("../data/combined_Robjects/SCE_QC_norm.rds")
ump <- read.csv("../data/combined_Robjects/UMAP_corrected.csv",stringsAsFactors=FALSE)[,-1]
sce <- sce[,ump$barcode] # this ensure the right order and cells for the umap graph below

# Graph
umapgraph <- readRDS("../data/combined_Robjects/UMAP_graphs.rds")
igr.cor <- umapgraph[["Corrected"]]

# Corrected PCA
# m.cor <- readRDS("../data/combined_Robjects/CorrectedPCA.rds")
# m.cor <- m.cor[colnames(sce),] # I should fix this to have to avoid doing this

# ---- Based on UMAP Graph ----
# Clustering with walktrap on the UMAP graph
set.seed(42)
cl.ump.wktrp.3 <- cluster_walktrap(igr.cor, steps=3)
set.seed(42)
cl.ump.wktrp.4 <- cluster_walktrap(igr.cor, steps=4)
set.seed(42)
cl.ump.wktrp.5 <- cluster_walktrap(igr.cor, steps=5)

# ---- Save ----
out <- data.frame("barcode"=colnames(sce),
		  "Walktrap3"=paste0("C",cl.ump.wktrp.3$membership),
		  "Walktrap4"=paste0("C",cl.ump.wktrp.4$membership),
		  "Walktrap5"=paste0("C",cl.ump.wktrp.5$membership))
write.csv(out,"../data/combined_Robjects/Compare_Clusters.csv")
