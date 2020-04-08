# Load Data
library(scran)
library(Matrix)
library(batchelor)
library(BiocParallel)
library(BiocSingular)
library(umap)
library(igraph)
source("../functions.R")


# Read in Tumorigenesis Data
sce <- readRDS("../../data/combined_Robjects/SCE_final.rds")
colnames(sce) <- sce$barcode # not sure why it wasn't like this in the first place
		      
# Read in Pregnancy Data
sce.g <- readRDS("../../data/combined_Robjects/SCE_final_Pregnancy.rds")

print("Read in data")
# Subset to same genes
genes <- intersect(rownames(sce),rownames(sce.g))
sce <- sce[genes,]
sce.g <- sce.g[genes,]

# Sub Sample for testing
# sce <- sce[,sample(ncol(sce),5000)]
# sce.g <- sce.g[,sample(ncol(sce.g),5000)]

# Scale normalization factors
mBatch <- batchelor::multiBatchNorm(sce,sce.g)
sce <- mBatch[[1]]
sce.g <- mBatch[[2]]
rm(mBatch)


# Extract ColData and combine
pD1 <- data.frame(colData(sce))
pD2 <- data.frame(colData(sce.g))
pD2$CellTypes <- paste0("G_",pD2$CellTypes)
pD2$Colors <- NA
pD2$Batch <- as.numeric(as.character(pD2$Batch)) + 3 
pD1$Batch <- as.numeric(as.character(pD1$Batch))

cols <- intersect(colnames(pD1),colnames(pD2))
cols <- cols[!grepl("UMAP",cols)] # One of them had the UMAP in there
pDout <- rbind(pD1[,cols],pD2[,cols])

# Split into the separate batches
sce1 <- sce[,sce$Batch==1]
sce2 <- sce[,sce$Batch==2]
sce3 <- sce[,sce$Batch==3]
sce4 <- sce.g[,sce.g$Batch==1]
sce5 <- sce.g[,sce.g$Batch==2]
# rm(sce)
rm(sce.g)
print("Prepared data")

#  Compute HVGs
# Compute highly variable genes per batch
dec.var1 <- modelGeneVar(sce1)
dec.var2 <- modelGeneVar(sce2)
dec.var3 <- modelGeneVar(sce3)
dec.var4 <- modelGeneVar(sce4)
dec.var5 <- modelGeneVar(sce5)


# Combine Var
combVar <- combineVar(dec.var1,dec.var2,dec.var3,dec.var4,dec.var5)
combVar <- combVar[rowData(sce)$KeepForHvg,]
hvgs <- getTopHVGs(combVar)
print("Combined Var")
# mnnCorrect
# param <- MulticoreParam(workers=4)

# Define merge order so that the pregnancy is mapped onto tumorigenesis
# The reasoning is that the tumorigenesis is the most heterogeneous, e.g. tumor cells
# If this is not the first one, these might (and will) be wrongly associated to what is the closest 
# In the other dataset
ord <- list(list(1:3),list(4:5))

print("Starting MNN")
set.seed(300)
mnncor <- batchelor::fastMNN(sce1,sce2,sce3,sce4,sce5,
			     #                      BPPARAM=param, # commented this out for singularity
		     k=20,
		     d=50,
		     merge.order=ord,
		     BSPARAM=IrlbaParam(deferred=TRUE),
		     #                      BNPARAM=AnnoyParam(),
		     #cos.norm=TRUE, # prob necesarry as I didn't use multibatchnorm
		     subset.row=hvgs,
		     correct.all=TRUE)
print("Done MNN")

m <- cbind(logcounts(sce1),logcounts(sce2),logcounts(sce3),
	   logcounts(sce4),logcounts(sce5))
rownames(pDout) <- pDout$barcode
pDout <- DataFrame(pDout[colnames(mnncor),])
colData(mnncor) <- pDout
assays(mnncor) <- list("logcounts"=m,
		       "reconstructed"=assays(mnncor)[["reconstructed"]])

# PCA
m.cor <- reducedDim(mnncor)

# Save correction
saveRDS(m.cor,"../../data/combined_Robjects/Pregnancy_CorrectedPCA.rds")
saveRDS(mnncor,"../../data/combined_Robjects/Pregnancy_FullMNN.rds")

# Compute UMAP for whole dataset
ump.cor <- umap(m.cor, random_state=42)
out <- data.frame("barcode"=rownames(m.cor))
out$UMAP1 <- ump.cor$layout[,1]
out$UMAP2 <- ump.cor$layout[,2]

igr.cor <- get_umap_graph(ump.cor)


# Save UMAP
saveRDS(igr.cor,"../../data/combined_Robjects/PregnancyUMAP_graph.rds")
write.csv(out,file="../../data/combined_Robjects/PregnancyUMAP_corrected.csv")
