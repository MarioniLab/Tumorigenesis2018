# Identifying doublets with doubletCluster
library(ggplot2)
library(cowplot)
library(wesanderson)
library(viridis)
theme_set(theme_cowplot())

# Read in data
out <- readRDS("backup.rds")

# Visualizing the clustering
cluster_plots <- lapply(names(out), function(nm) {
		   x <- out[[nm]]
		   p <- ggplot(x, aes(x=SubUMAP1, y=SubUMAP2, color=factor(Cluster))) +
			       geom_point() +
			       scale_color_manual(values=wes_palette("FantasticFox1",length(unique(x$Cluster)),type="continuous")) +
			       theme(legend.position="none", panel.background = element_rect(colour = "black", linetype = 1, size = 0.5),
				     axis.line = element_blank(),
				     axis.ticks = element_blank(),
				     axis.text = element_blank(),
				     axis.title = element_blank(),) +
			       ggtitle(paste0(nm," (n=",nrow(x),", k=",length(unique(x$Cluster)),")")) 
			   return(p)
})
plot_grid(plotlist=cluster_plots)

# Visualizing the doublet scores
score_plots <- lapply(names(out), function(nm) {
		   x <- out[[nm]]
		   p <- ggplot(x, aes(x=SubUMAP1, y=SubUMAP2, color=log2(DbltScore+1))) +
			       geom_point() +
			       scale_color_viridis() +
			       theme(legend.position="none", panel.background = element_rect(colour = "black", linetype = 1, size = 0.5),
				     axis.line = element_blank(),
				     axis.ticks = element_blank(),
				     axis.text = element_blank(),
				     axis.title = element_blank(),) +
			       ggtitle(paste0(nm," (n=",nrow(x),", k=",length(unique(x$Cluster)),")")) 
			   return(p)
})
plot_grid(plotlist=score_plots)


# Identifying doublets based on doubletCell
# Add sample index to data frames
outmod <- lapply(names(out), function(x) {
		dfs <- out[[x]]
		dfs$Sample <- x
		return(dfs)
})

outdf <- do.call(rbind,outmod)

# Compute median scores per cluster and sample
medscores <- aggregate(outdf$DbltScore, list(outdf$Sample, outdf$Cluster), median)
names(medscores) = c("sample", "cluster", "median.score")

mad_upper <- function(x){
  x <- x-median(x)
  return(mad(x[x>0], center = 0))
}

# Compute P-Values per sample based on normal null with mean of the median and sd of the mad
tests <- lapply(unique(outdf$Sample), function(x){
  sbset <- medscores[medscores$sample == x,]
  scores_sample <- outdf$DbltScore[outdf$Sample == x]
  sbset$p.value <- pnorm(sbset$median.score, mean = median(sbset$median.score), sd = mad_upper(sbset$median.score), lower.tail = FALSE)
  sbset$grp.FDR <- p.adjust(sbset$p.value,method="fdr") 
  return(sbset)
})

# Add a few stats
medscores <- do.call(rbind, tests)
medscores$fdr <- p.adjust(medscores$p.value, method = "fdr")
medscores$n.cells = sapply(1:nrow(medscores), function(row){
			       sum(outdf$Cluster == medscores$cluster[row] & outdf$Sample == medscores$sample[row])
})

medscores$frac.cells = sapply(1:nrow(medscores), function(row){
				  medscores$n.cells[row]/sum(medscores$n.cells[medscores$sample == medscores$sample[row]])
})

# Plot from Johnny with the significant clusters highlighted
ggplot(medscores, aes(x = sample, y = log2(median.score + 1), col = factor(cluster))) +
  geom_point(data = medscores[medscores$grp.FDR < 0.1,], size = 2) +
  geom_jitter(data = medscores[medscores$grp.FDR >= 0.1,], col = "darkgrey", size = 0.4, width = 0.2, height = 0) +
  theme(legend.position = "none") +
  #   scale_x_continuous(breaks = 1:max(medscores$sample)) +
  labs(x = "Sample", y = "log2(cluster median score + 1)")

#### Old Code
# plss <- list()
# for (i in seq_along(out)) {
#     test <- out[[i]]
#     m.s <- m.raw[,test$barcode]
#     dbls <- doubletCluster(m.s,test$Cluster)
#     crit1 <- dbls$N < 50
#     crit2 <- dbls$lib.size1 < 1 & dbls$lib.size2 < 1 & dbls$prop < 0.07
#     dblCluster <- rownames(dbls)[crit1 & crit2]

    # 
    #     ggplot(test, aes(x=factor(Cluster), y=log2(DbltScore+1))) +
    #         geom_violin() +
    #         geom_point() 

#     if (length(dblCluster)!=0) {
#     plss[[i]] <- ggplot(test, aes(x=SubUMAP1, y=SubUMAP2, color=Cluster %in% c(10))) +
#         geom_point() +
#         ggtitle(paste0("Sample ", names(out)[i], " has ", length(dblCluster), " Doublet Cluster")) 
# }
# }
