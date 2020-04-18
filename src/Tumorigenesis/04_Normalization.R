library(scran)
library(BiocSingular)
library(BiocParallel)
sce <- readRDS(file="../data/Robjects/SCE_QC.rds")

# Split into two batches
sce1 <- sce[,sce$Batch==1]
sce2 <- sce[,sce$Batch==2]
sce3 <- sce[,sce$Batch==3]

# Batch 1
set.seed(42)
clusters <- quickCluster(sce1, method="igraph", use.ranks=TRUE, d=50, BSPARAM=IrlbaParam(),BPPARAM=MulticoreParam(4),min.mean=0.01)
sce1 <- computeSumFactors(sce1, clusters=clusters)

plot(log10(sizeFactors(sce1)),log10(colSums(counts(sce1))),pch=19,xlab="Log(SizeFactors)",ylab="Log(LibrarySize)")

# Batch 2
set.seed(42)
clusters <- quickCluster(sce2, method="igraph", use.ranks=TRUE, d=50, BSPARAM=IrlbaParam(),BPPARAM=MulticoreParam(4),min.mean=0.01)
sce2 <- computeSumFactors(sce2, clusters=clusters)

plot(log10(sizeFactors(sce2)),log10(colSums(counts(sce2))),pch=19,xlab="Log(SizeFactors)",ylab="Log(LibrarySize)")

# Batch 3
set.seed(42)
clusters <- quickCluster(sce3, method="igraph", use.ranks=TRUE, d=50, BSPARAM=IrlbaParam(),BPPARAM=MulticoreParam(4),min.mean=0.01)
sce3 <- computeSumFactors(sce3, clusters=clusters)

plot(log10(sizeFactors(sce3)),log10(colSums(counts(sce3))),pch=19,xlab="Log(SizeFactors)",ylab="Log(LibrarySize)")

#multiBatchNorm
scemnorm <- batchelor::multiBatchNorm(sce1,sce2,sce3)
sceout <- cbind(scemnorm[[1]],scemnorm[[2]],scemnorm[[3]])
saveRDS(sceout,file="../data/Robjects/SCE_QC_norm.rds")
