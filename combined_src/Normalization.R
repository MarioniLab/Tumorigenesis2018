library(scran)
library(BiocSingular)
library(BiocParallel)
dataList <- readRDS(file="../data/combined_Robjects/ExpressionList_QC.rds")
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]

clusters <- quickCluster(m, method="igraph", use.ranks=TRUE, d=50, BSPARAM=IrlbaParam(),BPPARAM=MulticoreParam(4))
sce <- computeSumFactors(sce, clusters=clusters)
plot(log10(sizeFactors(sce)),log10(pD$UmiSums),pch=19,xlab="Log(SizeFactors)",ylab="Log(LibrarySize)")

sce <- normalize(sce)
m <- logcounts(sce)

out <- list()
out[["counts"]] <- m
out[["phenoData"]] <- pD
out[["featureData"]] <- fD
saveRDS(out,file="../data/combined_Robjects/ExpressionList_QC_norm.rds")
