# Load Data
library(scran)
library(Matrix)
library(dplyr)
library(BiocParallel)
library(BiocSingular)
library(BiocNeighbors)
library(igraph)
library(umap)

# Read in Data
sce <- readRDS("../data/combined_Robjects/SCE_QC_norm.rds")
		      
smps <- unique(sce$SampleID)
out <- bplapply(smps, function(smp) {
    # Subset to sample
    cells <- sce$SampleID==smp

    sce.sub <- sce[,cells]

    # Doublet Scores
    set.seed(42)
    scrs <- doubletCells(sce.sub,BSPARAM=IrlbaParam())

    # Quick SNN Graph on raw counts (!) as I want to cluster doublets together
    set.seed(42)
    igr <- buildSNNGraph(sce.sub,
			 BSPARAM=IrlbaParam(),
			 assay.type="counts",
			 BNPARAM=AnnoyParam()) # If slow maybe set a threshold on mean

    # Clustering on graph w high resolution
    cl <- cluster_walktrap(igr,steps=2) # steps reduced to 2 to get overly resolved clusters
    sce.sub$Cluster <- cl$membership

    # Highly variable genes for UMAP
    dec.var <- modelGeneVar(sce.sub)
    dec.var <- dec.var[rowData(sce.sub)$KeepForHvg,]
    hvg <- getTopHVGs(dec.var, prop=0.1) # prop 0.1 for speed purposes

    #UMAP and save
    ump <- umap(as.matrix(t(logcounts(sce.sub)[hvg,])), random_state=42)
    sce.sub$SubUMAP1 <- ump$layout[,1]
    sce.sub$SubUMAP2 <- ump$layout[,2]
    sce.sub$DbltScore <- scrs
    print(paste0("Done with: ",smp))
    out <- data.frame(colData(sce.sub))
    out <- out[,c("barcode","Cluster","SubUMAP1","SubUMAP2","DbltScore")]
    return(out)
},  BPPARAM=MulticoreParam(4))
names(out) <- smps
saveRDS(out,"../data/combined_Robjects/DoubletData.rds")
