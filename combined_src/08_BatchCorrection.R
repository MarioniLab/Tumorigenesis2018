# Load Data
library(scran)
library(ggplot2)
library(Matrix)
library(cowplot)
library(dplyr)
library(batchelor)
library(BiocParallel)
library(BiocSingular)
source("functions.R")

# Read in Data
dataList <- readRDS("../data/combined_Robjects/ExpressionList_QC_norm.rds")
# dataList <- subSample(dataList,cell.number=30000)
		      
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
rm(dataList)

# Remove Doublets 
qcMets <- read.csv("../data/combined_Robjects/QC_Part2.csv",row.names=1,stringsAsFactors=FALSE)
rmCells <- qcMets$barcode[qcMets$isDoubletFinal | qcMets$isRbc]

pD <- pD[!pD$barcode %in% rmCells,]

# Split into batch1 and batch2 

b1 <- pD[pD$Batch==1,"barcode"]
b2 <- pD[pD$Batch==2,"barcode"]

m1 <- m[,colnames(m) %in% b1]
m2 <- m[,colnames(m) %in% b2]

#  Compute HVGs
# Compute highly variable genes per batch and get overlap
hvg1 <- getHighVar(m1,get.var.out=TRUE)
hvg2 <- getHighVar(m2,get.var.out=TRUE)

# Remove genes that I don't consider for HVG incl Ribosomal genes and Mitochondrial genes
rmGenes <- fD$id[!fD$KeepForHvg]
combVar <- combineVar(hvg1,hvg2)
combVar$id <- rownames(combVar) 
combVar <- combVar[!(combVar$id %in% rmGenes),]

# take top 1000 for mnnCorrect
combVar <- arrange(data.frame(combVar), desc(bio))
hvg <- combVar$id[1:2000]

# mnnCorrect
param <- MulticoreParam(workers=2)
set.seed(300)

mnncor <- batchelor::fastMNN(m2,m1,
		     BPPARAM=param,
		     k=20,
		     BSPARAM=IrlbaParam(),
		     cos.norm=FALSE, # because I did multiBatchNorm shouldn't be necessary
		     subset.row=hvg)

# After correction
out.final <- reducedDim(mnncor)
out.final <- out.final[pD$barcode,] # reorder
saveRDS(out.final,"../data/combined_Robjects/CorrectedPCA.rds")
saveRDS(mnncor,"../data/combined_Robjects/FullMNN.rds")
