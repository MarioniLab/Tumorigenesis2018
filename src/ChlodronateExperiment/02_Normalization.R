library(scran)
library(BiocSingular)
library(BiocParallel)
sce <- readRDS(file="../../data/Tumorigenesis/Robjects/Chlodronate_SCE_QC.rds")

# Batch 1
set.seed(42)
clusters <- quickCluster(sce, method="igraph", use.ranks=TRUE, d=50, BSPARAM=IrlbaParam(),BPPARAM=MulticoreParam(4),min.mean=0.01)
sce <- computeSumFactors(sce, clusters=clusters)

# plot(log10(sizeFactors(sce1)),log10(colSums(counts(sce1))),pch=19,xlab="Log(SizeFactors)",ylab="Log(LibrarySize)")
#multiBatchNorm
saveRDS(sce,file="../../data/Tumorigenesis/Robjects/Chlodronate_SCE_QC_norm.rds")
