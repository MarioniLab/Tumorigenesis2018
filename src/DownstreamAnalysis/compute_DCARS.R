library(scran)
library(scater)
library(Matrix)
library(dplyr)
source("../TwoSampleDCARS.R")

sce.full <- readRDS("../../data/Integrated//Robjects/SCE_combined_final.rds")
sce.full$SubCellTypes <- factor(sce.full$CellTypes)
sce.full$CellTypes <- factor(sce.full$CellTypesFinal)
sce.full$Colors <- factor(sce.full$Colors)
sce.full$Batch <- as.character(sce.full$Batch)
sce.full$Groups <- factor(sce.full$Groups)
fkingp53_andtm <- c("WKBR76.5c","WKBR76.5f")

ptime <- read.csv("../../data/Integrated/Robjects/TumorTime.csv")[,-1]

pD <- data.frame(colData(sce.full))
pD <- dplyr::left_join(pD,ptime)
pD <- droplevels(pD)
colData(sce.full) <- DataFrame(pD)

sce.full <- sce.full[,sce.full$MajorGroups=="Epithelial" & !(sce.full$Condition %in% fkingp53_andtm)]

sce.full$Condition[grepl("WKBR",sce.full$Condition)] <- sce.full$ptimeBin[grepl("WKBR",sce.full$Condition)]
sce.full$Condition[grepl("CTRL",sce.full$Condition)] <- "WTYoung"
sce.full$Condition[grepl("CBLA1|CBLT",sce.full$Condition)] <- "WTOld"
sce.full$Condition <- factor(sce.full$Condition, levels=c("WTYoung","WTOld","4.5dG","9.5dG","14.5dG",
							  "1","2","3","4","5"))

sce <- sce.full[,sce.full$Experiment=="Tumorigenesis"]
rD <- rowData(readRDS("../../data/Tumorigenesis/Robjects/SCE.rds"))
rmgenes <- rownames(rD)[!rD$KeepForHvg]
rmgenes <- rownames(sce)[rownames(sce) %in% rmgenes]

# sce.lp <- sce.full[,sce.full$CellTypes %in% c("Lp","Avd") & sce.full$Condition!="WTOld"]
# cd <- data.frame(colData(sce.lp))
# cd <- droplevels(cd)
# colData(sce.lp) <- DataFrame(cd)
# dcars <- TwoSampleDCARS(sce.lp,
#                         assay="logcounts",
#                         PPI="Csn2",
#                         twoSample="Experiment",
#                         twoSampleRefLevel="Pregnancy",
#                         split=NULL,
#                         cor_method="spearman",
#                         minNonZeroSamples=10,
#                         niter=1000,
#                         ncores=4,
#                         parallel=TRUE)
# dcars <- saveRDS(dcars,"../../data/Downstream/DCARS_OUT.rds")

# The same for Lp Only
sce.lp <- sce.full[,sce.full$CellTypes %in% c("Lp") & sce.full$Condition!="WTOld"]
cd <- data.frame(colData(sce.lp))
cd <- droplevels(cd)
colData(sce.lp) <- DataFrame(cd)
dcars <- TwoSampleDCARS(sce.lp,
			assay="logcounts",
			PPI="Csn2",
			twoSample="Experiment",
			twoSampleRefLevel="Pregnancy",
			split=NULL,
			cor_method="spearman",
			minNonZeroSamples=10,
			niter=1000,
			ncores=4,
			parallel=TRUE)
saveRDS(dcars,"../../data/Downstream/DCARS_OUT_LpOnly.rds")

 sce.lp <- sce.full[,sce.full$CellTypes %in% c("Lp","Avd") & sce.full$Condition!="WTOld"]
 cd <- data.frame(colData(sce.lp))
 cd <- droplevels(cd)
 colData(sce.lp) <- DataFrame(cd)
 dcars <- TwoSampleDCARS(sce.lp,
                         assay="logcounts",
                         PPI="Cebpb",
                         twoSample="Experiment",
                         twoSampleRefLevel="Pregnancy",
                         split=NULL,
                        cor_method="spearman",
                         minNonZeroSamples=10,
                         niter=1000,
                         ncores=4,
                         parallel=TRUE)
 saveRDS(dcars,"../../data/Downstream/DCARS_OUT_Cebpb.rds")
