# Compute UMAP for whole dataste
library(dplyr)
library(ggplot2)
library(Matrix)
library(umap)
library(scran)
library(igraph)
library(BiocSingular)
source("../functions.R")

# ---- After Correction -----
m.cor <- readRDS("../../data/combined_Robjects/Pregnancy_CorrectedPCA.rds")
ump.cor <- umap(m.cor, random_state=42)
out <- data.frame("barcode"=rownames(m.cor))
out$UMAP1 <- ump.cor$layout[,1]
out$UMAP2 <- ump.cor$layout[,2]

igr.cor <- get_umap_graph(ump.cor)

saveRDS(igr.cor,"../../data/combined_Robjects/PregnancyUMAP_graph.rds")
write.csv(out,file="../../data/combined_Robjects/PregnancyUMAP_corrected.csv")
