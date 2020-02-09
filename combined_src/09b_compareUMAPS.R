# Compute UMAP for whole dataste
library(dplyr)
library(ggplot2)
library(Matrix)
library(umap)
library(scran)
library(igraph)
source("functions.R")

sce <- readRDS("../data/combined_Robjects/SCE_QC_norm.rds")

# Remove Doublets 
qcMets <- read.csv("../data/combined_Robjects/QC_Part2.csv",row.names=1,stringsAsFactors=FALSE)
rmCells <- qcMets$barcode[qcMets$isDoubletFinal | qcMets$isRbc]

sce <- sce[,!sce$barcode %in% rmCells]
# ---- After Correction -----
m.cor <- readRDS("../data/combined_Robjects/CorrectedPCA.rds")
m.cor <- m.cor[sce$barcode,]
min.dists <- c(0.1,0.25,0.5)
n.neighbors <- c(5,15,25,50,200)

out <- bplapply(n.neighbors, function(n.neighbor){
		    tmp <- lapply(min.dists, function(min.dist) {
				      ump <- umap(m.cor, random_state=42,
						  min_dist=min.dist,
						  n_neighbors=n.neighbor)
				      rtrn <- data.frame("barcode"=rownames(m.cor),
							 "UMAP1"=ump$layout[,1],
							 "UMAP2"=ump$layout[,2])
				      return(rtrn)
			    })
		    names(tmp) <- paste0("n_",n.neighbor,"_md_",min.dists)
		    return(tmp)}, BPPARAM=MulticoreParam(min(6L,length(n.neighbors))))
saveRDS(out,"../data/combined_Robjects/UMAP_Comparison.rds")
