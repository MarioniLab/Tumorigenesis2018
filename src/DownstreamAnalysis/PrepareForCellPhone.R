library(scran)
library(scater)
library(ggplot2)
library(cowplot)
library(Matrix)
library(umap)
library(dplyr)
library(viridis)
source("../functions.R")
theme_set(theme_cowplot())

sce.full <- readRDS("../../data/Integrated//Robjects/SCE_combined_final.rds")
sce.full$SubCellTypes <- factor(sce.full$CellTypes)
sce.full$CellTypes <- factor(sce.full$CellTypesFinal)
sce.full$Colors <- factor(sce.full$Colors)
sce.full$Batch <- as.character(sce.full$Batch)
sce.full$Groups <- factor(sce.full$Groups)
fkingp53_andtm <- c("WKBR76.5c","WKBR76.5f","WKBR61.4dTM","WKBR75.4bTM")
sce.full <- sce.full[, !(sce.full$Condition %in% fkingp53_andtm)]
ptime <- read.csv("../../data/Integrated/Robjects/TumorTime.csv")[,-1]

colData(sce.full) <- DataFrame(dplyr::left_join(data.frame(colData(sce.full)),ptime))
sce.full$Condition[grepl("WKBR",sce.full$Condition)] <- sce.full$ptimeBin[grepl("WKBR",sce.full$Condition)]
sce.full$Condition[grepl("CTRL",sce.full$Condition)] <- "WTYoung"
sce.full$Condition[grepl("CBLA1|CBLT",sce.full$Condition)] <- "WTOld"
sce.full$Condition <- factor(sce.full$Condition, levels=c("WTYoung","WTOld","4.5dG","9.5dG","14.5dG",
							  "1","2","3","4"))

rD <- rowData(readRDS("../../data/Tumorigenesis/Robjects/SCE.rds"))
outdr <- "../../data/Downstream/CellPhoneDB_out/"

## CellphoneDB

### Pregnancy
sce.p <- sce.full[,sce.full$Experiment=="Pregnancy"]
# Get count matrix with orthologs

library("biomaRt")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
 
mapping <- getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = rD$ID , mart = mouse, attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=T)
colnames(mapping)[1] <- "ID"
rD <- dplyr::left_join(data.frame(rD),mapping) 
rD <- rD[!is.na(rD$Gene.stable.ID.1),]
rD <- rD[rD$uniq %in% rownames(sce.full),]

## Prepare for CellPhoneDB
sce.out <- sce.p[rD$uniq,]
rownames(sce.out) <- rD$Gene.stable.ID.1
# SubSample For test
# set.seed(42)
# sce.out <- sce.out[,sample(ncol(sce.out),2000)]
Conds <- unique(sce.out$Condition)
for (Cond in Conds) {
    meta <- data.frame(colData(sce.out[,sce.out$Condition==Cond]))
    meta <- meta[,c("barcode","CellTypesFinal")]
    colnames(meta) <- c("Cell","cell_type")
    out <- as.matrix(logcounts(sce.out[,sce.out$Condition==Cond]))
    meta$Cell <- gsub("-",".",meta$Cell)
    colnames(out) <- meta$Cell
    write.table(out,paste0(outdr,"cellphonedb_count_",Cond,".txt"), sep="\t", quote=FALSE)
    write.table(meta,paste0(outdr,"cellphonedb_meta_",Cond,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}

### Tumor 
sce <- sce.full[,!is.na(sce.full$ptimeBin)]
# Get count matrix with orthologs
library("biomaRt")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
 
mapping <- getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = rD$ID , mart = mouse, attributesL = c("ensembl_gene_id"), martL = human, uniqueRows=T)
colnames(mapping)[1] <- "ID"
rD <- dplyr::left_join(data.frame(rD),mapping) 
rD <- rD[!is.na(rD$Gene.stable.ID.1),]
rD <- rD[rD$uniq %in% rownames(sce),]

## Prepare for CellPhoneDB
sce.out <- sce[rD$uniq,]
rownames(sce.out) <- rD$Gene.stable.ID.1
# SubSample For test
# set.seed(42)
# sce.out <- sce.out[,sample(ncol(sce.out),2000)]
bins <- c(1:4)
for (bin in bins) {
    meta <- data.frame(colData(sce.out[,sce.out$ptimeBin==bin]))
    meta <- meta[,c("barcode","CellTypesFinal")]
    colnames(meta) <- c("Cell","cell_type")
    out <- as.matrix(logcounts(sce.out[,sce.out$ptimeBin==bin]))
    meta$Cell <- gsub("-",".",meta$Cell)
    colnames(out) <- meta$Cell
    write.table(out,paste0(outdr,"cellphonedb_count_Bin",bin,".txt"), sep="\t", quote=FALSE)
    write.table(meta,paste0(outdr,"cellphonedb_meta_Bin",bin,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
