library(scran)
library(BiocSingular)
library(BiocParallel)
source("functions.R")
dataList <- readRDS(file="../data/combined_Robjects/ExpressionList_QC.rds")
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]

# Split into two batches
b1 <- pD[pD$Batch==1,"barcode"]
b2 <- pD[pD$Batch==2,"barcode"]

# Batch 1
m1 <- m[,colnames(m) %in% b1]
sce1 <- SingleCellExperiment(assays=list(counts=m1))
set.seed(42)
clusters <- quickCluster(m1, method="igraph", use.ranks=TRUE, d=50, BSPARAM=IrlbaParam(),BPPARAM=MulticoreParam(4))
sce1 <- computeSumFactors(sce1, clusters=clusters)

plot(log10(sizeFactors(sce1)),log10(colSums(m1)),pch=19,xlab="Log(SizeFactors)",ylab="Log(LibrarySize)")

# Batch 2
m2 <- m[,colnames(m) %in% b2]
sce2 <- SingleCellExperiment(assays=list(counts=m2))
set.seed(42)
clusters <- quickCluster(m2, method="igraph", use.ranks=TRUE, d=50, BSPARAM=IrlbaParam(),BPPARAM=MulticoreParam(4))
sce2 <- computeSumFactors(sce2, clusters=clusters)
plot(log10(sizeFactors(sce2)),log10(colSums(m2)),pch=19,xlab="Log(SizeFactors)",ylab="Log(LibrarySize)")


#multiBatchNorm
sce <- batchelor::multiBatchNorm(sce1,sce2)

# sce <- normalize(sce)
sf1 <- sizeFactors(sce[[1]])
sf2 <- sizeFactors(sce[[2]])
m1.norm <- logcounts(sce[[1]])
m2.norm <- logcounts(sce[[2]])

sfs <- c(sf1,sf2)

m.norm <- cbind(m1.norm,m2.norm)
names(sfs) <- colnames(m.norm)
m.norm <- m.norm[,pD$barcode]
sfs <- sfs[colnames(m.norm)]
add <- data.frame("barcode"=names(sfs),
		  "SizeFactor"=unname(sfs))

pD <- dplyr::left_join(pD,add)
out <- list()
out[["counts"]] <- m.norm
out[["phenoData"]] <- pD
out[["featureData"]] <- fD
saveRDS(out,file="../data/combined_Robjects/ExpressionList_QC_norm.rds")
