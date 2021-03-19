# Compute UMAP for whole dataste
library(dplyr)
library(ggplot2)
library(Matrix)
library(umap)
library(scran)
library(igraph)
library(BiocSingular)
source("../functions.R")

sce <- readRDS("../../data/Tumorigenesis/Robjects/SCE_QC_norm.rds")
# sce <- sce[,sample(ncol(sce),5000)]
# Remove Anything that failed second round of QC 
qcMets <- read.csv("../../data/Tumorigenesis/Robjects/QC_Part2.csv",row.names=1,stringsAsFactors=FALSE)
rmCells <- qcMets$barcode[qcMets$isDoubletFinal | qcMets$isRbc]

sce <- sce[,!colnames(sce) %in% rmCells]

var.dec <- modelGeneVar(sce)
var.dec <- var.dec[rowData(sce)$KeepForHvg,]
hvgs <- getTopHVGs(var.dec)
out <- data.frame("barcode"=colnames(sce))

set.seed(100)
pcs <- runPCA(t(logcounts(sce)[hvgs,]),BSPARAM=IrlbaParam(),rank=50)

# UMAP
pca.uncor <- pcs$x
rownames(pca.uncor) <- colnames(sce)
ump <- umap(pca.uncor, random_state=42)
out$UMAP1Uncor <- ump$layout[,1]
out$UMAP2Uncor <- ump$layout[,2]

igr.uncor <- get_umap_graph(ump)

# ---- After Correction -----
m.cor <- readRDS("../../data/Tumorigenesis/Robjects/CorrectedPCA.rds")
m.cor <- m.cor[colnames(sce),]
ump.cor <- umap(m.cor, random_state=42)
out$UMAP1 <- ump.cor$layout[,1]
out$UMAP2 <- ump.cor$layout[,2]

igr.cor <- get_umap_graph(ump.cor)


out2 <- list("Corrected"=igr.cor,
	    "Uncorrected"=igr.uncor)
saveRDS(out2,"../../data/Tumorigenesis/Robjects/UMAP_graphs.rds")
write.csv(out,file="../../data/Tumorigenesis/Robjects/UMAP_corrected.csv")
