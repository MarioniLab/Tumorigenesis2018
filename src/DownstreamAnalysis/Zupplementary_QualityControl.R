library(scran)
library(scater)
library(ggplot2)
library(cowplot)
library(dplyr)
library(edgeR)
library(Matrix)
library(ggrepel)
library(ggrastr)
library(viridis)
theme_set(theme_cowplot())
source("../functions.R")

sce <- readRDS("../../data/Integrated//Robjects/SCE_combined_final.rds")

sce$SampleID <- plyr::mapvalues(sce$SampleID,
				c("CBLA13","CBLT","CTRL_1","CTRL_2","CTRL_3"),
				c("WTOld_1","WTOld_2","WTYoung_1","WTYoung_2","WTYoung_3"),
				)

# Number of cells
pD <- data.frame(colData(sce))
smry <- group_by(pD, Experiment, Condition, SampleID) %>%
    summarize(nCells=n())

ncells <- ggplot(smry, aes(x=SampleID, y=nCells, fill=Condition)) +
    geom_bar(stat="identity",color="black") +
    facet_grid(~Condition, scales="free_x",space="free") +
    theme(panel.grid.major=element_line(colour="grey80",size=0.1,linetype="dashed"),
	  axis.text=element_text(size=14),
	  strip.text.x = element_text(size=12, face="bold") ,
	  strip.text.y = element_text(size=12, face="bold"),
	  legend.position="none",
	  strip.background = element_rect(colour="black", fill="white"),
	  axis.text.x=element_blank(),
	  axis.ticks.x=element_blank()) +
    scale_fill_brewer(palette="Paired") +
    ylab("Number of Cells") +
    xlab("")
   

umis <- ggplot(pD, aes(x=SampleID, y=UmiSums, fill=Condition)) +
    geom_violin(draw_quantiles=.5) +
    scale_y_log10() +
    facet_grid(~Condition, scales="free_x",space="free") +
    theme(panel.grid.major=element_line(colour="grey80",size=0.1,linetype="dashed"),
	  axis.text=element_text(size=14),
	  strip.text.x = element_blank() ,
	  strip.text.y = element_blank(),
	  strip.background = element_rect(colour="black", fill="white"),
	  legend.position="none",
	  axis.text.x=element_blank(),
	  axis.ticks.x=element_blank()) +
    scale_fill_brewer(palette="Paired") +
    ylab("# Total UMIs") +
    xlab("")

gds <- ggplot(pD, aes(x=SampleID, y=GenesDetected, fill=Condition)) +
    geom_violin(draw_quantiles=.5) +
    scale_y_log10() +
    facet_grid(~Condition, scales="free_x",space="free") +
    theme(panel.grid.major=element_line(colour="grey80",size=0.1,linetype="dashed"),
	  axis.text=element_text(size=14),
	  strip.text.x = element_blank() ,
	  strip.text.y = element_blank(),
#	  axis.text.x=element_text(size=14,angle=90),
	  axis.text.x=element_blank(),
	  axis.ticks.x=element_blank(),
	  strip.background = element_rect(colour="black", fill="white"),
	  legend.position="none") +
    scale_fill_brewer(palette="Paired") +
    ylab("# Genes detected")


# Celltype Composition
sumry <- group_by(pD, SampleID, Groups) %>%
    summarize(n=n()) %>%
    mutate(frac=n/sum(n)) %>%
    ungroup() 

comb <- table(sumry$Groups,sumry$SampleID) 
missing.comb <- which(comb==0, arr.ind=TRUE)
if (nrow(missing.comb)>0) {
    add.missing <- data.frame("Groups"=rownames(missing.comb),
		      "SampleID"=colnames(comb)[missing.comb[,"col"]],
		      "frac"=0,
		      "n"=0)
    sumry <- rbind(sumry,add.missing)
}

sumry$Groups <- factor(sumry$Groups, levels=c("TumorEpithel","Basal","Luminal",
					      "Myeloid","Lymphoid","Fibroblast",
					      "Stroma"))

sumry$Colors <- plyr::mapvalues(as.character(sumry$Groups),
				as.character(unique(pD$Groups)),
				as.character(unique(pD$GroupColors)))

sumry <- arrange(sumry,Groups)
add <- unique(pD[,c("SampleID","Condition")])
sumry <- dplyr::left_join(sumry,add)
sumry$Condition

barplt <- ggplot(sumry, aes(x=SampleID, y=frac*100, fill=Groups)) +
    geom_bar(stat="identity",color="black") +
    scale_fill_manual(values=unique(sumry$Colors)) +
    facet_grid(~Condition, scales="free_x",space="free") +
    theme(panel.grid.major=element_line(colour="grey80",size=0.1,linetype="dashed"),
	  axis.text=element_text(size=14),
	  strip.text.x = element_blank() ,
	  strip.text.y = element_blank(),
	  axis.text.x=element_text(size=14,angle=90),
	  strip.background = element_rect(colour="black", fill="white"), 
	  legend.position="none") +
    xlab("Condition") +
    ylab("Percentage of Cells") 

plot_grid(ncells,umis,gds,barplt,nrow=4,rel_heights=c(1,1,1,2))

# Marker Genes
smrzd <- aggregateAcrossCells(sce,ids=sce$CellTypesFinal,
			      use_exprs_values="logcounts",
			      average=TRUE)

ordr <- c("Hs","Hsp","Lp","Avd","LpBsl","Bsl1","Bsl2","BslG","Tm",
	  "Fb1","Fb2","Fb3","Fb4","Fb5","Fb6","Fb7","Fb8","Fb9","Cafs",
	  "Ec1","Ec2","Ec3","Ec4","Pericytes1","Pericytes2","Swc",
	  "cDC1","cDC2","migDC","pDC","MdC3","MdC1","MdC2","MastCells",
	  "Neutrophils","Mo3","Mo1","Mo2","Tam1","Tam2",
	  "CD81","CD82","CD83","CTLs","CD42","CD41","Tregs",
	  "CyclingT","NK","ILCs","BCells","PlasmaCells")
genes <- c("Epcam","Lgals7","Sfn",
	   "Krt8",
	   "Pgr","Prlr","Esr1","Cited1",
	   "Aldh1a3","Fcgbp","Lurap1l",
	   "Elf5","Csn2","Wap",
	   "Lgals7",
	   "Trp63","Gjb4",
	   "Cntn2","Krt5","Oxtr","Nrtn","Arc","Acta2",#LpBsl
	   "Col4a1","Nrg1",
	   "Oas2","Asz1",
	   "Col3a1","Dcn","Dpep1","Dpt","Lum","Pdgfra",
	   "Peg3","Lamc3","Zim1",
	   "Creb5",
	   "Ifi205",
	   "Dpt","C3",
	   "Mfap5","Pi16",
	   "Gdf10","Adh1",
	   "Smoc2","Apod","Gpc3",
	   "Sfrp5","Plscr1","Itgb4","C2",
	   "Tac1",
	   "Srgap1","Nptx2",
	   "Tnc","Col12a1","Col8a1",
	   "Fabp4","Apold1",
	   "Emcn","Pecam1",
	   "Nts","Lyve1",
	   "Vwf","Selp","Aqp1",
	   "Gkn3","Sema3g",
	   "Gm13889",
	   "Notch3",
	   "Mpz","Mbp",
	   "Ptprc","Cd52","Laptm5",
	   "Tyrobp","H2-Aa",
	   "Xcr1","Wdfy4",
	   "Tnip3","Cd209a","Ccl22","Il4i1",
	   "Ccr7",
	   "Siglech",
	   "Ear2",
	   "Chil3",
	   "Adgre4",
	   "Cpa3",
	   "S100a8",
	   "Cd14","Adgre1","Mmp12",
	   "Mrc1",
	   "Cxcl2",
	   "Arg1","Spp1","Ms4a7",
	   "Hcst","Cd2",
	   "Cd3d",
	   "Cd8a","Sell",
	   "Dapl1","Limd2","Cd27",
	   "Xcl1","Gzmb",
	   "Cd4","Icos","Tmem64","Frat2","Satb1","Foxp3","Ctla4",
	   "Mki67",
	   "Eomes","Klra8","Klra4","Klri2",
	   "Rora","Il7r","Cxcr6","Cd163l1","Il23r",
	   "Cd19","Pax5","H2-Ob","Fcmr","Jchain","Mzb1","Derl3")
mtx <- logcounts(smrzd)
mtx <- mtx/rowMax(mtx)
mtx <- mtx[genes,ordr]
library(viridis)
library(pheatmap)
colnames(mtx) <- renameForPlot(colnames(mtx))
pheatmap(t(mtx),
	 color=inferno(n=100,begin=0,end=1),
	 cluster_rows=FALSE,
	 cluster_cols=FALSE,
	 gaps_row=c(9,19,26,40),
	 fontsize=8
	 )
