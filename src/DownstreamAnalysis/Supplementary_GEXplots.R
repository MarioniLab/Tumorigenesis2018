library(scran)
library(scater)
library(cowplot)
library(dplyr)
library(edgeR)
library(Matrix)
library(ggrepel)
library(ggrastr)
library(viridis)
library(pheatmap)
require(reshape2)
theme_set(theme_cowplot())
source("../functions.R")

sce <- readRDS("../../data/Integrated//Robjects/SCE_combined_final.rds")
pD <- data.frame(colData(sce))
pD$CellTypesFinal <- renameForPlot(pD$CellTypesFinal)


#genes <- c("Lifr","Lif","Fgf1","Il33")
gene <- "Mki67"
pD$Expr <- logcounts(sce)[gene,]

ggplot(pD, aes(x=CellTypesFinal, y=Expr)) +
    geom_violin(scale="width", draw_quantiles=.5) +
    facet_grid(Condition~Groups,scales="free_x",space="free") +
    theme(axis.text.x=element_text(angle=45,hjust=1)) +
    ggtitle(gene)

# Plot for reviewer's question
######
sce.epi <- sce[,sce$MajorGroups=="Epithelial"]
pD.epi <- data.frame(colData(sce.epi))
pD.epi$CellTypesFinal <- renameForPlot(pD.epi$CellTypesFinal)

pD.epi <- pD.epi[,c("barcode","CellTypesFinal")]
genes <- c("Sox10","Cebpb","Elf5","Foxc1")
#genes <- c("Aldh1a3")

for (gene in genes) {
    pD.epi[,gene] <- logcounts(sce.epi)[gene,]
}

pD.epi <- melt(pD.epi)

ggplot(pD.epi, aes(x=CellTypesFinal, y=value)) +
    geom_violin(scale="width", draw_quantiles=.5) +
    facet_wrap(~variable,scales="free_y") +
    theme_pub() +
    theme(axis.text.x=element_text(angle=45,hjust=1)) +
    theme(strip.background = element_rect(colour="white", fill="white")) +
    ylab("Expression") +
    xlab("")
##########
ggsave("../../data/figures/Supplementary/TFs.pdf")

pD <- group_by(pD,Condition,Groups,CellTypesFinal) %>%
    summarize(nPos=sum(Expr>0),
	      numb=n()) %>%
    mutate(nFrac=nPos/numb)

ggplot(pD, aes(x=CellTypesFinal,y=nFrac)) +
    geom_bar(stat="identity") +
    facet_grid(Groups~Condition,scales="free_x")

sce.sub <- sce[,sce$Condition=="2"]
bubblePlot(logcounts(sce.sub),markers=c("Mki67","Top2a"),sce.sub$CellTypesFinal)
library(schex)

sce.sub <- make_hexbin(sce.sub, nbins = 100, 
		   dimension_reduction = "UMAP", use_dims=c(1,2))
pltGEX <- function(sce, gene) {
    plot_hexbin_gene(sce, gene=gene, type="logcounts",
		     action="mean") +
	theme_void() +
	theme(legend.position="none") +
	scale_fill_viridis(option="C") +
	ggtitle("")
}
pltGEX(sce.sub,"Mki67")


sce.sub <- sce[]
smrzd <- aggregateAcrossCells(sce,
			      id=colData(sce)[,c("Condition")],
			      average=TRUE,
			      use_exprs_values="logcounts")
ordr <- arrange(data.frame(colData(smrzd)),Condition)$Condition
genes <- c("Lif","Fgf1","Spp1","Tnfsf4","Tnfrsf4","Arg1")
mtx <- logcounts(smrzd)
mtx <- mtx[genes,ordr[-10]]
mtx <- mtx/rowMax(mtx)
pheatmap(mtx,
	 color=inferno(n=100,begin=0,end=1),
	 cluster_rows=FALSE,
	 cluster_cols=FALSE,
#	 gaps_col=c(2),
	 fontsize=8
	 )

sce.sub <- sce[,sce$CellTypesFinal=="Avd" & sce$Condition %in% c("4.5dG","9.5dG","14.5dG","1","2","3","4")]
gene <- "Spp1"
pD <- droplevels(data.frame(colData(sce.sub)))
pD$Expr <- logcounts(sce.sub)[gene,]
ggplot(pD, aes(x=Condition, y=Expr)) +
    geom_violin(scale="width", draw_quantiles=.5) +
    theme(axis.text.x=element_text(angle=45,hjust=1)) +
    ylab("Spp1") 


#sce.sub <- sce[]
smrzd <- aggregateAcrossCells(sce.sub,
			      id=droplevels(colData(sce.sub)[,c("Condition")]),
			      average=TRUE,
			      use_exprs_values="logcounts")

ordr <- arrange(data.frame(colData(smrzd)),Condition)$Condition
genes <- c("Ccl5")
mtx <- logcounts(smrzd)
mtx <- mtx[genes,]
mtx <- mtx/max(mtx)
pheatmap(t(mtx),
	 color=inferno(n=100,begin=0,end=1),
	 cluster_rows=FALSE,
	 cluster_cols=FALSE,
#	 gaps_col=c(2),
	 fontsize=8
	 )

sce.sub <- sce[, sce$MajorGroups=="Epithelial" & sce$Condition %in% c("WTOld","WTYoung","4.5dG","9.5dG","14.5dG","1","2","3","4","5")]


#sce.sub <- sce[]
smrzd <- aggregateAcrossCells(sce.sub,
			      id=droplevels(colData(sce.sub)[,c("Condition","CellTypesFinal")]),
			      average=TRUE,
			      use_exprs_values="logcounts")

ordr <- arrange(data.frame(colData(smrzd)),Condition)$Condition

genes <- c("Spp1")
mtx <- logcounts(smrzd)
mtx <- mtx[genes,]
pD <- data.frame(colData(smrzd))
pD$Expression <- mtx
pD$CellTypesFinal <- renameForPlot(pD$CellTypesFinal)

ggplot(pD, aes(x=Condition, y=CellTypesFinal, fill=Expression)) +
    geom_tile() +
    scale_fill_viridis(option="B") +
    theme_pub() +
    theme(axis.title.y=element_blank(),
	panel.grid.major=element_blank(),
	axis.text.x = element_text(angle = 45, hjust = 1),
	legend.position="none") +
    xlab("")
ggsave("../../data/figures/Supplementary//Spp1Expression_wTumour.pdf",width=5.69,height=3.75)

pD <- pD[pD$CellTypesFinal!="Tm" & pD$Condition!="5",]
ggplot(pD, aes(x=Condition, y=CellTypesFinal, fill=Expression)) +
    geom_tile() +
    scale_fill_viridis(option="B") +
    theme_pub() +
    theme(axis.title.y=element_blank(),
	panel.grid.major=element_blank(),
	axis.text.x = element_text(angle = 45, hjust = 1),
	legend.position="none") +
    xlab("")

ggsave("../../data/figures/Figure4/Spp1Expression.pdf",width=5.69,height=3.75)


sce.sub <- sce[,sce$Condition %in% c("WTOld","WTYoung","4.5dG","9.5dG","14.5dG","1","2","3","4","5")]

#sce.sub <- sce[]
smrzd <- aggregateAcrossCells(sce.sub,
			      id=droplevels(colData(sce.sub)[,c("Condition","CellTypesFinal")]),
			      average=FALSE,
			      use_exprs_values="logcounts")

ordr <- arrange(data.frame(colData(smrzd)),Condition)$Condition
genes <- c("Mki67")
mtx <- logcounts(smrzd)
mtx <- mtx[genes,]
pD <- data.frame(colData(smrzd))
pD$Expression <- mtx
pD$CellTypesFinal <- renameForPlot(pD$CellTypesFinal)
pD <- pD[pD$Condition %in% c("WTOld","1","2","3","4","5"),]
ggplot(pD, aes(x=Condition, y=CellTypesFinal, fill=Expression)) +
    geom_tile() +
    scale_fill_viridis(option="B") +
    theme_pub() +
    theme(axis.title.y=element_blank(),
	panel.grid.major=element_blank(),
	axis.text.x = element_text(angle = 45, hjust = 1),
	axis.text.y = element_text(size=14),
	legend.position="none") +
    xlab("")
ggsave("../../data/figures/Supplementary/Ki67.pdf",width=5.67,height=9.5)
