# Compute UMAP for whole dataste
library(dplyr)
library(ggplot2)
library(Matrix)
library(umap)
library(scran)
library(igraph)

umps <- readRDS("../data/combined_Robjects/UMAP_Comparison.rds")
comb <- lapply(umps,function(x) do.call(cbind,x))
comb <- do.call(cbind,comb)
comb <- comb[,!grepl("barcode",colnames(comb))]
comb$barcode <- rownames(comb)
umps.lt <- reshape::melt(comb)
umps.lt$variable <- as.character(umps.lt$variable)
umps.lt$Neighbors <- unlist(lapply(strsplit(umps.lt$variable,"_"),function(x) x[2]))
umps.lt$MinDist <- unlist(lapply(strsplit(umps.lt$variable,"_"),function(x) x[4]))
umps.lt$MinDist <- unlist(lapply(strsplit(umps.lt$MinDist,"\\."),function(x) paste0(x[1],x[2])))
umps.lt$variable <- unlist(lapply(strsplit(umps.lt$variable,"\\."),function(x) x[3]))
fplot <- spread(umps.lt,variable,value)

fplot$Neighbors <- factor(fplot$Neighbors,levels=c("5","15","25","50","200"))
fplot <- fplot[fplot$Neighbors!="5",]
library(cowplot)
theme_set(theme_cowplot())
library(ggrastr)
ggplot(fplot, aes(x=UMAP1, y=UMAP2)) +
    geom_point_rast() +
    facet_grid(MinDist~Neighbors, scales="free")
ggsave("UMAP_Comparison.pdf",width=16,height=12)
