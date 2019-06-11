library(RColorBrewer)
library(reshape2)
library(plyr)
library(dplyr)
library(scran)
library(Rtsne)
library(ggplot2)
library(Matrix)
library(igraph)
library(dynamicTreeCut)
library(cluster)
library(cowplot)
library(scater)
library(irlba)
source("functions.R")


dataList <- readRDS("../data/Robjects/ExpressionList_QC.rds")
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
rm(dataList)
sfs <- read.csv("../data/Robjects/SizeFactors.csv")
cluster <- read.csv("../data/Robjects/CellTypes.csv",row.names=1)

pD <- inner_join(pD,cluster)
pD <- left_join(pD,sfs)

# Gene and cell filtering
m <- m[,pD$barcode]

# Norm and Log
m <- log2(t(t(m)/pD$sf)+1)
m <- as(m,"dgCMatrix")

compUmap <-  function(m, fD, frctHvg=10) {
	require(umap)
	keep <- rowMeans(m)>0.01
	m <- m[keep,]

	# Highly variable genes
	hvg <- getHighVar(m,get.var.out=TRUE, supress.plot=TRUE)
	hvg <- hvg[rownames(hvg) %in% fD$uniqnames[fD$KeepForHvg],]
	hvg <- rownames(hvg[order(hvg$bio, decreasing=TRUE),])[1:(nrow(hvg)/frctHvg)]

	set.seed(300)
	ump <- umap(as.matrix(t(m[hvg,])),min_dist=0.5,,metric="pearson")

	out <- data.frame("barcode"=colnames(m),
			  "UMAP1"=ump$layout[,1],
			  "UMAP2"=ump$layout[,2])
	return(out)
}


grps <- unique(pD$Groups)

out <- data.frame()
for (grp in grps) {
    m.sub <- m[,pD$barcode[pD$Groups %in% grp]]
    tmp <- compUmap(m.sub, fD)
    out <- rbind(out,tmp)
    print(paste("Finished",grp))
}

colnames(out)[2:3] <- paste("Group",colnames(out)[2:3],sep=".")
pD.out <- right_join(pD,out)
# ggplot(pD.out, aes(x=Group.UMAP1, y=Group.UMAP2, color=CellType)) +
#     geom_point() +
#     facet_wrap(~Groups) +
#     scale_color_manual(values=levels(pD.out$Color)) 

out <- data.frame()
celltps <- names(table(pD$CellType)[table(pD$CellType)>100])
for (celltp in celltps) {
    m.sub <- m[,pD$barcode[pD$CellType %in% celltp]]
    tmp <- compUmap(m.sub, fD, frctHvg=20)
    out <- rbind(out,tmp)
}
colnames(out)[2:3] <- paste("CellType",colnames(out)[2:3],sep=".")
pD.out <- left_join(pD.out,out)


res <- pD.out[,c("barcode","Group.UMAP1","Group.UMAP2","CellType.UMAP1","CellType.UMAP2")]
write.csv(file="../data/Robjects/Substructure.csv",res)
