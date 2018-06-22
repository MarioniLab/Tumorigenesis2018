# Estimate size factors using scran
library(dplyr)
library(scran)
library(Matrix)
library(Rtsne)
library(limSolve)
source("functions.R")

dataList <- readRDS("../data/Robjects/ExpressionList_QC.rds")
# set.seed(300)
# dataList <- subSample(dataList, cell.number=5000)
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
rm(dataList)

# Gene and cell filtering
m <- m[fD$keep,pD$PassAll]
pD <- pD[pD$PassAll,]
fD <- fD[fD$keep,]
rownames(m) <- as.vector(rownames(m))

clusters <- quickCluster(m,method="igraph")
minSize <- min(table(clusters))
param <- MulticoreParam(workers=4)
pD$sf <- computeSumFactors(m, sizes=seq(20,min(100,minSize),5),clusters=clusters,positive=TRUE, BPPARAM=param)

plot(log10(colSums(m))~log10(pD$sf),main="Library Size versus Size Factors")

out <- pD[,c("barcode","sf")]
write.csv(file="../data/Robjects/SizeFactors.csv",out,row.names=FALSE)
