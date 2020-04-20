library(scran)
library(scater)

# Combined SCE
sce <- readRDS("../../data/Integrated/Robjects/Pregnancy_FullMNN.rds")
# UMAP
ump <- read.csv("../../data/Integrated/Robjects/PregnancyUMAP_corrected.csv",stringsAsFactors=FALSE)[,-1]
# Read in TumorTime
ptime <- read.csv("../../data/Tumorigenesis/Robjects/TumorTime.csv",stringsAsFactors=FALSE)[,-1]
ptime$ptime <- round(ptime$ptime,1)
# Read in Clustering results
cluster <- read.csv("../../data/Integrated/Robjects/PregnancyClusters.csv",stringsAsFactors=FALSE)[,-1]

# Put everything in pD
pD <- data.frame(colData(sce))
rmcols <- c("PassLibSize", "PassGenesDetected", "PassViability", "PassTred",
	    "PassAll", "Cluster", "Colors", "Barcode")
pD <- pD[,! (colnames(pD) %in% rmcols)]
colnames(pD)[colnames(pD)=="CellTypes"] <- "PreCellTypes"
pD <- dplyr::left_join(pD,ump)
reducedDim(sce,type="UMAP") <- pD[,c("UMAP1","UMAP2")]
pD <- dplyr::left_join(pD,ptime)
pD <- dplyr::left_join(pD,cluster)

#Randomize order for plotting
set.seed(42)
sce <- sce[,sample(ncol(sce),ncol(sce))]
rownames(pD) <- pD$barcode
pD <- pD[sce$barcode,]

# Sort out meta-data
pD$SampleID <- as.character(pD$SampleID)
pD$Condition <- as.character(pD$Condition)
# For Group Column
ctrls <- grep("CB|CT",unique(pD$SampleID),value=TRUE)
gest <- grep("dG",unique(pD$SampleID),value=TRUE)
wkbr <- grep("WKBR",unique(pD$SampleID),value=TRUE)
mp <- c(ctrls,gest,wkbr)
names(mp) <- c(rep("Control",length(ctrls)),
	       rep("Gestation",length(gest)),
	       rep("WKBR_MG",length(wkbr)))

pD$Group <- plyr::mapvalues(pD$SampleID,mp,names(mp))
pD$Group[grepl("TM",pD$SampleID)] <- "WKBR_Tumor"

# For Condition Column
wkbr <- as.character(ptime$Condition)
names(wkbr) <- paste0("TT_",ptime$ptime)
ctrls <- grep("CB|CT",unique(pD$Condition),value=TRUE)
names(ctrls) <- rep("Control",length(ctrls))
mp <- c(ctrls,wkbr)

pD$Condition <- plyr::mapvalues(pD$Condition,mp,names(mp))
mxdsrt <- gtools::mixedsort(names(wkbr))
pD$Condition <- factor(pD$Condition, levels=c("Control","4.5dG",
					     "9.5dG","14.5dG",
					     mxdsrt))

stopifnot(identical(rownames(pD),colnames(sce)))
colData(sce) <- DataFrame(pD)
saveRDS(sce,file="../../data/Integrated/Robjects/SCE_combined.rds")
