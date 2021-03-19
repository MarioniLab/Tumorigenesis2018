# Compute UMAP for whole dataste
library(dplyr)
library(ggplot2)
library(Matrix)
library(umap)
library(scran)
library(igraph)
library(BiocSingular)
source("../functions.R")

sce <- readRDS("../../data/Tumorigenesis/Robjects/Chlodronate_SCE_QC_norm.rds")
sce <- scater::logNormCounts(sce)

var.dec <- modelGeneVar(sce)
var.dec <- var.dec[rowData(sce)$KeepForHvg,]
hvgs <- getTopHVGs(var.dec)
out <- data.frame("barcode"=colnames(sce))

set.seed(100)
pcs <- runPCA(t(logcounts(sce)[hvgs,]),BSPARAM=IrlbaParam(),rank=50)

# UMAP
pca <- pcs$x
rownames(pca) <- colnames(sce)
ump <- umap(pca, random_state=42)
out$UMAP1 <- ump$layout[,1]
out$UMAP2 <- ump$layout[,2]

#igr <- get_umap_graph(ump)

# out <- list("Corrected"=igr.cor,
#             "Uncorrected"=igr.uncor)
# saveRDS(out2,"../../data/Tumorigenesis/Robjects/UMAP_graphs.rds")
write.csv(out,file="../../data/Tumorigenesis/Robjects/Chlodronate_UMAP_corrected.csv")
