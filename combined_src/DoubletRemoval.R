# Load Data
library(scran)
library(ggplot2)
library(Matrix)
library(cowplot)
library(dplyr)
library(batchelor)
library(BiocParallel)
library(BiocSingular)
library(viridis)
library(igraph)
library(umap)
source("functions.R")

# Read in Data
dataList <- readRDS("../data/combined_Robjects/ExpressionList_QC_norm.rds")
# dataList <- subSample(dataList,cell.number=50000)
		      
m.raw <- readRDS("../data/combined_Robjects/ExpressionList_QC.rds")[["counts"]]
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
rownames(pD) <- pD$barcode
fD <- dataList[["featureData"]]
rm(dataList)

smps <- unique(pD$SampleID)
out <- bplapply(smps, function(smp) {
    # Subset to sample
    cells <- pD[pD$SampleID==smp,"barcode"]

    m.raw.sub <- m.raw[,cells]
    m.sub <- m[,cells]
    m.sub <- m.sub[rowMeans(m.sub)>0.01,]
    pD.sub <- pD[cells,]

    # Doublet Scores
    set.seed(2206)
    scrs <- doubletCells(m.raw.sub,BSPARAM=IrlbaParam(),size.factors.norm=pD.sub$SizeFactor)

    # Clustering on all genes
    igr <- buildSNNGraph(m.sub, BSPARAM=IrlbaParam())
    cl <- cluster_walktrap(igr,steps=2) # steps reduced to 2 to get overly resolved clusters
    pD.sub$Cluster <- cl$membership

    # Highly variable genes for UMAP
    hvg <- getHighVar(m.sub,get.var.out=TRUE, supress.plot=TRUE)
    hvg <- hvg[rownames(hvg) %in% fD$uniqnames[fD$KeepForHvg],]
    hvg <- rownames(hvg[order(hvg$bio, decreasing=TRUE),])[1:(nrow(hvg)/10)]

    ump <- umap(as.matrix(t(m.sub[hvg,])), min_dist=0.5, metric="pearson", random_state=42)
    pD.sub$SubUMAP1 <- ump$layout[,1]
    pD.sub$SubUMAP2 <- ump$layout[,2]
    pD.sub$DbltScore <- scrs
    print(paste0("Done with: ",smp))
    out <- pD.sub[,c("barcode","Cluster","SubUMAP1","SubUMAP2","DbltScore")]
    return(out)
},  BPPARAM=MulticoreParam())
saveRDS(out,"backup.rds")
