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

igr.uncor <- get_umap_graph(ump)

# ---- After Correction -----
m.cor <- readRDS("../data/combined_Robjects/CorrectedPCA.rds")
m.cor <- m.cor[pD$barcode,]
ump <- umap(m.cor, random_state=42)
pD$UMAP1 <- ump$layout[,1]
pD$UMAP2 <- ump$layout[,2]

igr.cor <- get_umap_graph(ump)


out <- list("Corrected"=igr.cor,
	    "Uncorrected"=igr.uncor)
saveRDS(out,"../data/combined_Robjects/UMAP_graphs.rds")
write.csv(pD[,c("barcode","UMAP1","UMAP2","UMAP1Uncor","UMAP2Uncor")],file="../data/combined_Robjects/UMAP_corrected.csv")
