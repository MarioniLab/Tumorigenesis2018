library(randomForest)
library(scran)
sce <- readRDS("../../data/Tumorigenesis/Robjects/Chlodronate_SCE_QC_norm.rds")
tree <- readRDS("../../data/Integrated/Robjects/Tree.rds")
ump <- read.csv("../../data/Tumorigenesis/Robjects/Chlodronate_UMAP_corrected.csv",stringsAsFactors=FALSE)[,-1]
sce <- scater::logNormCounts(sce)
pD <- data.frame(colData(sce))
pD <- dplyr::left_join(pD,ump)
pD$ChlodronateExperiment <- pD$Condition %in% c("WKBR96.1d","WKBR96.1e",
						"WKBR89.4h","WKBR89.4j",
						"WKBR89.4i")
pred <- predict(tree,t(logcounts(sce)))

pD$CellTypes <- pred
pD$CellTypes <- factor(pD$CellTypes)

X <- model.matrix(~0+pD$CellTypes)
colnames(X) <- levels(pD$CellTypes)
X <- t(t(X)/colSums(X))
celltyp_ctrs <- data.frame(crossprod(X, as.matrix(pD[,c("UMAP1","UMAP2")])))
celltyp_ctrs$CellType <- rownames(celltyp_ctrs)

library(ggplot2)
p0 <- ggplot(pD, aes(x=UMAP1, y=UMAP2)) +
    geom_point() +
    ggtitle("Cluster") +
    geom_label(data=celltyp_ctrs, aes(x=UMAP1, y=UMAP2, label=CellType, color=NULL)) +
    theme(legend.position="none") +
    facet_wrap(~ChlodronateExperiment)
p0

cluster.counts <- table(pD$CellTypes, pD$Condition)
cluster.counts <- cluster.counts[rowMeans(cluster.counts)>10,]

library(dplyr)
sumry <- group_by(pD, Batch, Condition, ChlodronateExperiment) %>%
    summarize()

sumry$ChlodronateExperiment <- as.character(sumry$ChlodronateExperiment)

# EdgeR
library(edgeR)
# Set up DGE
y.ab <- DGEList(cluster.counts)
y.ab$samples$Condition <- colnames(y.ab)
y.ab$samples <- left_join(y.ab$samples,sumry)
colnames(y.ab) <- y.ab$samples$Condition

# I am using calcNormFactors to correct for compositional biases introdcued by one cluster changing a lot (in this case Av/Lp)
# This assumes that the majority of clusters don't change, which is fair given it's a highly heterogeneous dataset (40 cls) 
# And it's still the same tissue after all and does not include the tumor samples
# y.ab <- calcNormFactors(y.ab)

des <- model.matrix(~ ChlodronateExperiment,y.ab$samples)

y.ab <- estimateDisp(y.ab, des, trend="none")
plotBCV(y.ab,cex=1)

fit.ab <- glmQLFit(y.ab, des, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
summary(fit.ab$df.prior)
plotQLDisp(fit.ab, cex=1)

res <- glmQLFTest(fit.ab, coef=ncol(fit.ab$design))
summary(decideTests(res,p.value=0.1))
tab <- topTags(res,n=40)$table
# tab <- tab[tab$FDR<0.1,]
# tab
# addCol <- data.frame("CellTypes"=levels(pD$CellTypes),
#                      "Colors"=levels(cold$Colors))
# rownames(addCol) <- addCol$CellTypes
# addCol <- addCol[rownames(tab),]
# tab$Colors <- addCol$Colors
tab$CellTypes <- rownames(tab)
# tab$Colors <- factor(tab$Colors, levels=tab$Colors)
tab$CellTypes <- factor(tab$CellTypes, levels=tab$CellTypes)
library(ggrepel)
ggplot(tab, aes(x=logFC, y=-log10(FDR),color=CellTypes)) +
    geom_point(size=4) +
    geom_hline(yintercept=1) +
    geom_label_repel(data=tab[tab$FDR<0.5,], label=rownames(tab[tab$FDR<0.5,]),aes(color=NULL)) +
    #     scale_color_manual(values=levels(tab$Colors)) +
    theme(legend.position="none")
