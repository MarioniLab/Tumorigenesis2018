# Load Data
library(scran)
library(Matrix)
library(batchelor)
library(BiocParallel)
library(BiocSingular)

# Read in Data
sce <- readRDS("../../data/Tumorigenesis/Robjects/SCE_QC_norm.rds")
		      
# Remove Anything that failed second round of QC 
qcMets <- read.csv("../../data/Tumorigenesis/Robjects/QC_Part2.csv",row.names=1,stringsAsFactors=FALSE)
rmCells <- qcMets$barcode[qcMets$isDoubletFinal | qcMets$isRbc]

sce <- sce[,!colnames(sce) %in% rmCells]
# Split into batch1 and batch2 

sce1 <- sce[,sce$Batch==1]
sce2 <- sce[,sce$Batch==2]
sce3 <- sce[,sce$Batch==3]

#  Compute HVGs
# Compute highly variable genes per batch
dec.var1 <- modelGeneVar(sce1)
dec.var2 <- modelGeneVar(sce2)
dec.var3 <- modelGeneVar(sce3)

# Combine Var
combVar <- combineVar(dec.var1,dec.var2,dec.var3)
stopifnot(identical(rownames(combVar),rownames(sce)))
combVar <- combVar[rowData(sce)$KeepForHvg,]
hvgs <- getTopHVGs(combVar)

# mnnCorrect
param <- MulticoreParam(workers=4)

set.seed(300)
mnncor <- batchelor::fastMNN(sce1,sce2,sce3,
		     BPPARAM=param,
		     k=20,
		     d=50,
		     BSPARAM=IrlbaParam(deferred=TRUE),
		     #                      BNPARAM=AnnoyParam(),
		     cos.norm=FALSE, # because I did multiBatchNorm shouldn't be necessary
		     subset.row=hvgs)

# After correction
out.final <- reducedDim(mnncor)
saveRDS(out.final,"../../data/Tumorigenesis/Robjects/CorrectedPCA.rds")
saveRDS(mnncor,"../../data/Tumorigenesis/Robjects/FullMNN.rds")
