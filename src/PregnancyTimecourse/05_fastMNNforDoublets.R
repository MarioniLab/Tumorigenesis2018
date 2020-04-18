# Load Data
library(scran)
library(Matrix)
library(batchelor)
library(umap)
library(BiocParallel)
library(BiocSingular)
library(BiocNeighbors)

# Read in Data
sce <- readRDS("../data/Robjects/SCE_QC_norm.rds")
# sce <- sce[,sample(1:ncol(sce),20000)]
# Split into batch1 and batch2 

sce1 <- sce[,sce$Batch==1]
sce2 <- sce[,sce$Batch==2]

#  Compute HVGs
# Compute highly variable genes per batch
dec.var1 <- modelGeneVar(sce1)
dec.var2 <- modelGeneVar(sce2)
# Combine Var
combVar <- combineVar(dec.var1,dec.var2)
stopifnot(identical(rownames(combVar),rownames(sce)))
combVar <- combVar[rowData(sce)$KeepForHvg,]
hvgs <- getTopHVGs(combVar)

# mnnCorrect
param <- MulticoreParam(workers=4)

set.seed(300)
mnncor <- batchelor::fastMNN(sce2,sce1,
		     BPPARAM=param,
		     k=20,
		     d=50,
		     BSPARAM=IrlbaParam(deferred=TRUE),
		     BNPARAM=AnnoyParam(),
		     cos.norm=FALSE, # because I did multiBatchNorm shouldn't be necessary
		     subset.row=hvgs)

# Extract PCA
pca.cor <- reducedDim(mnncor)

# Compute UMAP 
ump <- umap(pca.cor, random_state=42)
reducedDim(mnncor,type="umap") <- ump$layout[,1:2]

# Save
saveRDS(mnncor,"../data/Robjects/Corrected_preDoublet.rds")
