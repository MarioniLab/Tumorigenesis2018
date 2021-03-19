# This script is just to clean up the SCE object and add some color definitions
library(scran)
library(scater)

sce <- readRDS("../../data/Integrated/Robjects/SCE_combined_celltypes.rds")

# p53 WT animals
p53wt <- c("WKBR76.5c","WKBR76.5f")
sce <- sce[, !(sce$Condition %in% p53wt)]

# Add Tumor stages
ptime <- read.csv("../../data/Integrated/Robjects/TumorTime.csv")[,-1]
pD <- data.frame(colData(sce))
pD <- dplyr::left_join(pD,ptime)
pD <- droplevels(pD)

# Put back into SCE object
colData(sce) <- DataFrame(pD)

# Define Conditions
sce$Condition[grepl("WKBR",sce$Condition)] <- sce$ptimeBin[grepl("WKBR",sce$Condition)]
sce$Condition[grepl("CTRL",sce$Condition)] <- "WTYoung"
sce$Condition[grepl("CBLA1|CBLT",sce$Condition)] <- "WTOld"
sce$Condition <- factor(sce$Condition, levels=c("WTYoung","WTOld","4.5dG","9.5dG","14.5dG",
							  "1","2","3","4","5"))

sce$Groups <- factor(sce$Groups)
sce$MajorGroups <- factor(sce$MajorGroups)

# Batch is numeric, turn to character to prevent modelling as continuous
sce$Batch <- as.character(sce$Batch)

# Add Colors
sce$CellTypesFinal <- factor(sce$CellTypesFinal)

library(RColorBrewer)
# Eptihelium
epithelial <- unique(sce$CellTypesFinal[sce$MajorGroups=="Epithelial"])
getPalette <- colorRampPalette(brewer.pal(9, "Blues"))
col.epithelial <- getPalette(length(epithelial)+1)[2:(length(epithelial)+1)]
names(col.epithelial) <- epithelial
col.epithelial[names(col.epithelial)=="Tm"] <- "#800080"

# Myeloid
myeloid <- unique(sce$CellTypesFinal[sce$Groups=="Myeloid"])
getPalette <- colorRampPalette(brewer.pal(9, "Oranges"))
col.myeloid <- getPalette(length(myeloid)+1)[2:(length(myeloid)+1)]
names(col.myeloid) <- myeloid

# Lymphoid
lymphoid <- unique(sce$CellTypesFinal[sce$Groups=="Lymphoid"])
getPalette <- colorRampPalette(brewer.pal(9, "Reds"))
col.lymphoid <- getPalette(length(lymphoid)+1)[2:(length(lymphoid)+1)]
names(col.lymphoid) <- lymphoid

# Fibroblasts
fibroblasts <- unique(sce$CellTypesFinal[sce$MajorGroups=="Fibroblast"])
getPalette <- colorRampPalette(brewer.pal(9, "Greens"))
col.fibroblasts <- getPalette(length(fibroblasts)+1)[2:(length(fibroblasts)+1)]
names(col.fibroblasts) <- fibroblasts

# Stroma
stroma <- unique(sce$CellTypesFinal[sce$MajorGroups=="Stroma"])
col.stroma <- brewer.pal(n=length(stroma)+1,"Greys")[2:(length(stroma)+1)]
names(col.stroma) <- stroma

colors <- c(col.epithelial, col.myeloid, col.lymphoid, col.fibroblasts,
	    col.stroma)

sce$GroupColors <- plyr::mapvalues(sce$Groups,
				levels(sce$Groups),
				c("#00BFFF","#37A055","#3787C0","#F78757","#FF4E03",
				  "#737373","#800080"))

saveRDS(sce,"../../data/Integrated/Robjects/SCE_combined_final.rds")
