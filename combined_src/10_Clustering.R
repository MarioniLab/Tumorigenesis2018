# Load Data
dataList <- readRDS("../data/combined_Robjects/ExpressionList_QC_norm.rds")
pD <- dataList[["phenoData"]]
rownames(pD) <- pD$barcode
rm(dataList)
ump <- read.csv("../data/combined_Robjects/UMAP_corrected.csv",stringsAsFactors=FALSE)
pD <- pD[ump$barcode,] # this ensure the right order and cells for the umap graph below

library(igraph)
umapgraph <- readRDS("../data/combined_Robjects/UMAP_graphs.rds")
igr.cor <- umapgraph[["Corrected"]]
set.seed(42)
cl.cor <- cluster_walktrap(igr.cor,steps=3)
pD$Cluster <- as.factor(paste0("C",cl.cor$membership))
write.csv(file="../data/combined_Robjects/Clusters.csv",pD[,c("barcode","Cluster")])
