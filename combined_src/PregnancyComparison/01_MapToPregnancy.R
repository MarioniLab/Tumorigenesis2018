# Load Data
library(scran)
library(Matrix)
library(batchelor)
library(BiocParallel)
library(BiocSingular)

# Read in Tumorigenesis Data
sce <- readRDS("../../data/combined_Robjects/SCE_final.rds")
colnames(sce) <- sce$barcode # not sure why it wasn't like this in the first place
		      
# Read in Pregnancy Data
sce.g <- readRDS("../../../PregnancyTimecourse2019/data/Robjects/SCE_final.rds")

genes <- intersect(rownames(sce),rownames(sce.g))
sce <- sce[genes,]
sce.g <- sce.g[genes,]


pD1 <- data.frame(colData(sce))
pD2 <- data.frame(colData(sce.g))
pD2$CellTypes <- paste0("G_",pD2$CellTypes)
pD2$Colors <- NA
pD2$Batch <- as.numeric(as.character(pD2$Batch)) + 3 
pD1$Batch <- as.numeric(as.character(pD1$Batch))

cols <- intersect(colnames(pD1),colnames(pD2))
cols <- cols[!grepl("UMAP",cols)]
pDout <- rbind(pD1[,cols],pD2[,cols])

# Split into batch1 and batch2 
sce1 <- sce[,sce$Batch==1]
sce2 <- sce[,sce$Batch==2]
sce3 <- sce[,sce$Batch==3]
sce4 <- sce.g[,sce.g$Batch==1]
sce5 <- sce.g[,sce.g$Batch==2]
# rm(sce)
rm(sce.g)

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

# mnnCorrect
param <- MulticoreParam(workers=4)

set.seed(300)
mnncor <- batchelor::fastMNN(sce4,sce5,sce1,sce2,sce3,
		     BPPARAM=param,
		     k=20,
		     d=50,
		     BSPARAM=IrlbaParam(deferred=TRUE),
		     #                      BNPARAM=AnnoyParam(),
		     cos.norm=TRUE, # prob necesarry as I didn't use multibatchnorm
		     subset.row=hvgs,
		     correct.all=TRUE)

m <- cbind(logcounts(sce4),logcounts(sce5),logcounts(sce1),
	   logcounts(sce2),logcounts(sce3))
rownames(pDout) <- pDout$barcode
pDout <- DataFrame(pDout[colnames(mnncor),])
colData(mnncor) <- pDout
assays(mnncor) <- list("logcounts"=m,
		       "reconstructed"=assays(mnncor)[["reconstructed"]])

# After correction
out.final <- reducedDim(mnncor)
saveRDS(out.final,"../../data/combined_Robjects/Pregnancy_CorrectedPCA.rds")
saveRDS(mnncor,"../../data/combined_Robjects/Pregnancy_FullMNN.rds")
