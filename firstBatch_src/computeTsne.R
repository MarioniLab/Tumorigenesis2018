# Estimate size factors using scran
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(scran)
library(Rtsne)
library(ggplot2)
library(Matrix)
library(cowplot)
library(irlba)
source("functions.R")

dataList <- readRDS("../data/Robjects/ExpressionList_QC.rds")
pD.full <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
sfs <- read.csv("../data/Robjects/SizeFactors.csv")
pD.full <- right_join(pD.full, sfs)
conds <- list("38","41","46","53",
	      c("38","41"),c("38","41","53"),
	      c("38","41","46"),c("38","41","46","53"),
	      "nontumor")
batchCor <- c("raw")
plotlist <- list()
	      
for (bc in batchCor) {
    for (i in 1:length(conds)) {
	cond <- conds[[i]]

	if (bc=="raw") {
	m <- dataList[["counts"]]
	rownames(m) <- as.vector(rownames(m))
	} else {
	m <- dataList[["counts"]]
	}
	if (cond=="nontumor") {
	    cells <- pD.full$barcode[pD.full$PassAll & pD.full$Condition!="WKBR61.4dT"]
	} else {
	cells <- pD.full$barcode[pD.full$Age %in% cond & pD.full$PassAll] 
	}
	# Gene and cell filtering
	m <- m[,cells]
	keep <- rowMeans(m)>0.01 
	pD <- pD.full[pD.full$barcode %in% cells,]
	m <- m[keep,]

	#Norm and log
	m <- log2(t(t(m)/pD$sf)+1)

	# Highly variable genes
	hvg <- getHighVar(m,get.var.out=TRUE, supress.plot=TRUE)
	hvg <- hvg[rownames(hvg) %in% fD$id[fD$KeepForHvg],]
	hvg <- rownames(hvg[order(hvg$bio, decreasing=TRUE),])[1:(nrow(hvg)/10)]

	# Compute tSNE 
	fPCA <- m[hvg,]
	pca <- prcomp_irlba(t(fPCA), n=50, center=FALSE, scale.=FALSE)
	pca <- pca$x
	set.seed(300)
	tsn <- Rtsne(pca,perplexity=50,pca=FALSE)
	pD$tSNE1 <- tsn$Y[,1]
	pD$tSNE2 <- tsn$Y[,2]

	# plot
	pD <- arrange(pD, prcntMito)
	plotlist[[paste(cond,collapse=".")]] <- ggplot(pD, aes(x=tSNE1, y=tSNE2, color=Condition)) +
	    geom_point() +
	    theme_classic()

	# tSNE
	out <- pD[,c("barcode","tSNE1","tSNE2")]
	filename <- paste0("../data/Robjects/tSNE_",paste(cond,collapse="-"),"_",bc,".csv")
	write.csv(out,file=filename,row.names=FALSE,quote=FALSE)
    }
}

pdf("tSNEs_byAge.pdf",height=12,width=12)
plot_grid(plotlist=plotlist,nrow=2,ncol=2)
dev.off()
