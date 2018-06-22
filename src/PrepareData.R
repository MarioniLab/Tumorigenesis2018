# Script to prepare cellranger data for downstream analysis
library(plyr)
library(dplyr)
library(reshape2)

# ---- ReadData ----

# Read in output dropletUtils package
cDat <- readRDS("../data/Robjects/CountMatrix.rds")
pDat <- data.frame(barcode=colnames(cDat))
fDat <- read.table("../data/CellRangerData/SampleA12/outs/filtered_gene_bc_matrices/mm10/genes.tsv",stringsAsFactors=FALSE)
colnames(fDat) <- c("id","symbol")
rownames(fDat) <- fDat$id
fDat <- fDat[rownames(cDat),]

#reduce size of matrix
keep <- rowSums(cDat) > 1
cDat <- cDat[keep,]
fDat <- fDat[keep,]

# ---- Formatting ----

# conversion into sample names
conv <- read.csv("../data/miscData/SLX-15931.HV5MKBBXX.s_8.contents.csv", stringsAsFactors=FALSE)
pDat$barcode <- as.character(pDat$barcode)
conv$Sample.name <- gsub("-","",conv$Sample.name) 

# Add more info to phenotype Data
pDat <- mutate(pDat, SeqID=substr(barcode,18,nchar(barcode))) %>%
        mutate(SampleID=mapvalues(SeqID,conv$Barcode,
				  conv$Sample.name)) %>%
	mutate(Age=mapvalues(SampleID,conv$Sample.name,
			     conv$Age)) %>%
	mutate(Replicate=as.factor(substr(SampleID,nchar(SampleID),nchar(SampleID))),
	       Condition=as.factor(substr(SampleID,1,nchar(SampleID)-1)))
    
pDat$Replicate[pDat$Replicate=="M"] <- 1
pDat$Replicate[pDat$Replicate=="G"] <- 1

pDat$AnimalHadTumor <- pDat$Condition %in% c("WKBR61.4dT","WKBR61.4dM")
# Add more info to the feature Data
mitoGenes <- read.table("../data/miscData/MitoGenes.txt")
tfCheck <- read.table("../data/miscData/TFcheckpoint_WithENSID.tsv",
		header=TRUE, sep="\t")
# Ribosmal genes
library(biomaRt)
gos <- c("GO:0003735","G0:0005840","GO:0015935","GO:0015934") # constituents of ribosome
ensembl  <- useMart("ensembl",dataset="mmusculus_gene_ensembl") 
gene.data <- getBM(attributes=c('ensembl_gene_id', 'go_id'), mart=ensembl)
gene.data <- gene.data[gene.data$go_id %in% gos,]

fDat$Ribosomal <- fDat$id %in% gene.data$ensembl_gene_id

# Mitochondrial
fDat$Mitochondrial <- fDat$id %in% mitoGenes$V1

# TranscriptionFactor
fDat$TranscriptionFactor <- fDat$id %in% tfCheck$ensembl_gene_id

# Keep for hvg

fDat$KeepForHvg <- !(fDat$Ribosomal | fDat$Mitochondrial)


# Save data
stopifnot(identical(rownames(fDat),as.character(rownames(cDat))) & identical(colnames(cDat),pDat$barcode))
DataList <- list("phenoData"=pDat, "featureData"=fDat, "counts"=cDat)
saveRDS(DataList,file="../data/Robjects/ExpressionList.rds")
