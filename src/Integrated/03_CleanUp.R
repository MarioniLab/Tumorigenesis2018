library(scran)
library(scater)

# Combined SCE
sce <- readRDS("../../data/Integrated/Robjects/Pregnancy_FullMNN.rds")
# UMAP
ump <- read.csv("../../data/Integrated/Robjects/PregnancyUMAP_corrected.csv",stringsAsFactors=FALSE)[,-1]

cluster <- read.csv("../../data/Integrated/Robjects/PregnancyClusters.csv",stringsAsFactors=FALSE)[,-1]
# Put everything in pD
pD <- data.frame(colData(sce))
rmcols <- c("PassLibSize", "PassGenesDetected", "PassViability", "PassTrend",
	    "PassAll", "Barcode","Cluster","CellTypes","Groups","MajorGroups")
pD <- pD[,! (colnames(pD) %in% rmcols)]
pD <- dplyr::left_join(pD,ump)
pD <- dplyr::left_join(pD,cluster)
reducedDim(sce,type="UMAP") <- pD[,c("UMAP1","UMAP2")]

#Randomize order for plotting
set.seed(42)
sce <- sce[,sample(ncol(sce),ncol(sce))]
rownames(pD) <- pD$barcode
pD <- pD[sce$barcode,]

# Sort out meta-data
pD$SampleID <- as.character(pD$SampleID)
pD$Condition <- as.character(pD$Condition)

#For Sample Column
#Some WKBR samples were split on two lanes, hence two SeqIDs per animal
#The SampleID will reflect animals, so this will be rm
pD$SampleID[grepl("WKBR",pD$SampleID)] <- substr(pD$SampleID[grepl("WKBR",pD$SampleID)],1,(nchar(pD$SampleID[grepl("WKBR",pD$SampleID)])-1))
#Same for CBLT
pD$SampleID[grepl("CBLT",pD$SampleID)] <- substr(pD$SampleID[grepl("CBLT",pD$SampleID)],1,(nchar(pD$SampleID[grepl("CBLT",pD$SampleID)])-1))

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

stopifnot(identical(rownames(pD),colnames(sce)))
colData(sce) <- DataFrame(pD)
saveRDS(sce,file="../../data/Integrated/Robjects/SCE_combined.rds")
