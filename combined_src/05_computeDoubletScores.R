# Load Data
library(scran)
library(Matrix)
library(dplyr)
library(BiocParallel)
library(BiocSingular)
library(BiocNeighbors)
library(umap)

# Read in Data
sce <- readRDS("../data/combined_Robjects/SCE_QC_norm.rds")
		      
smps <- unique(sce$SampleID)
out <- bplapply(smps, function(smp) {
    # Subset to sample
    cells <- sce$SampleID==smp

    sce.sub <- sce[,cells]

    # Doublet Scores
    # Highly variable genes
    dec.var <- modelGeneVar(sce.sub)
    dec.var <- dec.var[rowData(sce.sub)$KeepForHvg,]
    hvgs <- getTopHVGs(dec.var, prop=0.5)  # this is for the sake of speed

    set.seed(42)
    scrs <- doubletCells(sce.sub,BSPARAM=IrlbaParam(),
			 subset.row=hvgs)

    # Quick SNN Graph 
    set.seed(42)
    igr <- buildSNNGraph(sce.sub,
			 BSPARAM=IrlbaParam(),
			 assay.type="logcounts", # for my own memory as I previously did it on counts
			 subset.row=hvgs,
			 k=10,
			 BNPARAM=AnnoyParam()) # If slow maybe set a threshold on mean
    # Quick KNNs
    set.seed(42)
    pc.approx <- runPCA(t(counts(sce.sub)), rank=50, BSPARAM=IrlbaParam())
    knn <- findKNN(pc.approx$x, k=10, BNPARAM=AnnoyParam())

    # Clustering on graph
    cl <- igraph::cluster_walktrap(igr, steps=2)
    sce.sub$Cluster <- cl$membership

    #UMAP and save
    ump <- umap(as.matrix(t(logcounts(sce.sub)[hvgs,])), random_state=42)
    sce.sub$SubUMAP1 <- ump$layout[,1]
    sce.sub$SubUMAP2 <- ump$layout[,2]
    sce.sub$DbltScore <- scrs
    print(paste0("Done with: ",smp))
    out <- colData(sce.sub)
    out$KNNs <- knn$index
    out <- out[,c("barcode","Cluster","SubUMAP1","SubUMAP2","DbltScore","KNNs")]
    return(out)
},  BPPARAM=MulticoreParam(4))
names(out) <- smps
saveRDS(out,"../data/combined_Robjects/DoubletData.rds")
