# Load Data
library(scran)
library(ggplot2)
library(Matrix)
library(cowplot)
library(dplyr)
library(batchelor)
library(BiocParallel)
library(BiocSingular)
library(viridis)
library(igraph)
library(umap)
source("functions.R")

# Read in Data
dataList <- readRDS("../data/combined_Robjects/ExpressionList_QC_norm.rds")
# dataList <- subSample(dataList,cell.number=30000)
		      
m.raw <- readRDS("../data/combined_Robjects/ExpressionList_QC.rds")[["counts"]]
m <- dataList[["counts"]]
pD <- dataList[["phenoData"]]
fD <- dataList[["featureData"]]
rm(dataList)

out <- list()
plts <- list()
smps <- unique(pD$SampleID)
for (smp in smps) {
    # Subset to sample
    cells <- pD[pD$SampleID==smp,"barcode"]

    m.raw.sub <- m.raw[,colnames(m) %in% cells]
    m.sub <- m[,colnames(m) %in% cells]
    m.sub <- m.sub[rowMeans(m.sub)>0.01,]
    pD.sub <- pD[pD$barcode %in% cells,]

    # Doublet Scores
    set.seed(2206)
    scrs <- doubletCells(m.raw.sub,BSPARAM=IrlbaParam(),size.factors.norm=pD.sub$SizeFactor)

    # Clustering on all genes
    igr <- buildSNNGraph(m.sub, BSPARAM=IrlbaParam(), BPPARAM=MulticoreParam(4))
    cl <- cluster_walktrap(igr)
    pD.sub$Cluster <- cl$membership

    # Highly variable genes for UMAP
    hvg <- getHighVar(m.sub,get.var.out=TRUE, supress.plot=TRUE)
    hvg <- hvg[rownames(hvg) %in% fD$uniqnames[fD$KeepForHvg],]
    hvg <- rownames(hvg[order(hvg$bio, decreasing=TRUE),])[1:(nrow(hvg)/10)]
    ump <- umap(as.matrix(t(m.sub[hvg,])), min_dist=0.5, metric="pearson", random_state=42)
    pD.sub$SubUMAP1 <- ump$layout[,1]
    pD.sub$SubUMAP2 <- ump$layout[,2]
    pD.sub$DbltScore <- scrs
    pD.sub <- arrange(pD.sub, DbltScore)
    plts[[smp]] <- ggplot(pD.sub, aes(x=UMAP1, y=UMAP2,color=DbltScore)) +
		    geom_point() +
		    scale_color_viridis() +
		    ggtitle(smp)
    out[[smp]] <- pD.sub[,c("barcode","Cluster","SubUMAP1","SubUMAP2","DbltScore")]
    print(paste0("Done with: ",smp))
}
saveRDS(out,"backup.rds")

# Identifying doublets with doubletCluster
plss <- list()
for (i in seq_along(out)) {
    test <- out[[i]]
    m.s <- m.raw[,test$barcode]
    dbls <- doubletCluster(m.s,test$Cluster)
    crit1 <- dbls$N < 50
    crit2 <- dbls$lib.size1 < 1 & dbls$lib.size2 < 1 & dbls$prop < 0.07
    dblCluster <- rownames(dbls)[crit1 & crit2]

    # 
    #     ggplot(test, aes(x=factor(Cluster), y=log2(DbltScore+1))) +
    #         geom_violin() +
    #         geom_point() 

    if (length(dblCluster)!=0) {
    plss[[i]] <- ggplot(test, aes(x=SubUMAP1, y=SubUMAP2, color=Cluster %in% c(10))) +
	geom_point() +
	ggtitle(paste0("Sample ", names(out)[i], " has ", length(dblCluster), " Doublet Cluster")) 
}
}

plot_grid(plotlist=plss[1:8])
plot_grid(plotlist=plss[9:17])



# Identifying doublets based on doubletCell
outmod <- lapply(names(out), function(x) {
		dfs <- out[[x]]
		dfs$Sample <- x
		return(dfs)
})

outdf <- do.call(rbind,outmod)

medscores <- aggregate(outdf$DbltScore, list(outdf$Sample, outdf$Cluster), median)
names(medscores) = c("sample", "cluster", "median.score")

mad_upper <- function(x){
  x <- x-median(x)
  return(mad(x[x>0], center = 0))
}

tests <- lapply(unique(outdf$Sample), function(x){
  sbset <- medscores[medscores$sample == x,]
  scores_sample <- outdf$DbltScore[outdf$Sample == x]
  sbset$p.value = pnorm(sbset$median.score, mean = median(sbset$median.score), sd = mad_upper(sbset$median.score), lower.tail = FALSE)
  return(sbset)
})

medscores <- do.call(rbind, tests)
medscores$fdr = p.adjust(medscores$p.value, method = "fdr")
medscores$n.cells = sapply(1:nrow(medscores), function(row){
  sum(clusters == medscores$cluster[row] & outdf$Sample == medscores$sample[row])
})

medscores$frac.cells = sapply(1:nrow(medscores), function(row){
  medscores$n.cells[row]/sum(medscores$n.cells[medscores$sample == medscores$sample[row]])
})



ggplot(medscores, aes(x = sample, y = log2(median.score + 1), col = factor(cluster))) +
  geom_point(data = medscores[medscores$fdr < 0.1,], size = 2) +
  geom_jitter(data = medscores[medscores$fdr >= 0.1,], col = "darkgrey", size = 0.4, width = 0.2, height = 0) +
  theme(legend.position = "none") +
  #   scale_x_continuous(breaks = 1:max(medscores$sample)) +
  labs(x = "Sample", y = "log2(cluster median score + 1)")
