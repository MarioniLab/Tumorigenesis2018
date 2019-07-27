# Script to prepare cellranger data for downstream analysis
library(plyr)
library(dplyr)
library(reshape2)

# ---- ReadData ----

# Read in output dropletUtils package
cDat <- readRDS("../data/combined_Robjects/CountMatrix.rds")
pDat <- data.frame(barcode=colnames(cDat))
fDat <- read.table("../data/CellRangerOutput/genes.tsv",stringsAsFactors=FALSE)
colnames(fDat) <- c("id","symbol")
rownames(fDat) <- fDat$id
fDat <- fDat[rownames(cDat),]

#reduce size of matrix
keep <- rowSums(cDat) > 1
cDat <- cDat[keep,]
fDat <- fDat[keep,]

# ---- Formatting ----

# conversion into sample names
conv.b1 <- read.csv("../data/misc/SLX-15931.HV5MKBBXX.s_8.contents.csv", stringsAsFactors=FALSE)
conv.b2 <- read.csv("../data/misc/SLX-18124.HCJKGDRXX.s_1.contents.csv", stringsAsFactors=FALSE)
conv <- rbind(conv.b1,conv.b2)
pDat$barcode <- as.character(pDat$barcode)
conv$Sample.name <- gsub("-","",conv$Sample.name) 

# Add more info to phenotype Data
pDat <- mutate(pDat, SeqID=substr(barcode,18,nchar(barcode))) %>%
        mutate(SampleID=mapvalues(SeqID,conv$Barcode,
				  conv$Sample.name)) %>%
	mutate(Age=mapvalues(SampleID,conv$Sample.name,
			     conv$Age)) %>%
	mutate(Chip=as.factor(substr(SampleID,nchar(SampleID),nchar(SampleID))),
	       Animal=substr(SampleID,1,nchar(SampleID)-1),
	       Condition=as.factor(substr(SampleID,1,nchar(SampleID)-1)))
    
pDat$Animal <- gsub("TM|MG","",pDat$Animal)
pDat$Batch <- ifelse(pDat$Chip=="3",2,1)

pDat$AnimalHadTumor <- pDat$Animal %in% c("WKBR61.4d","WKBR75.4b")
# Add more info to the feature Data
mitoGenes <- read.table("../data/misc/MitoGenes.txt")
tfCheck <- read.table("../data/misc/TFcheckpoint_WithENSID.tsv",
		header=TRUE, sep="\t")
# Ribosmal genes
library(biomaRt)
gos <- c("GO:0003735","G0:0005840","GO:0015935","GO:0015934") # constituents of ribosome
ensembl  <- useMart("ensembl",dataset="mmusculus_gene_ensembl",host="useast.ensembl.org") # at the time of writing the ensembl.org had issues
gene.data <- getBM(attributes=c('ensembl_gene_id', 'go_id'), mart=ensembl)
gene.data <- gene.data[gene.data$go_id %in% gos,]

fDat$Ribosomal <- fDat$id %in% gene.data$ensembl_gene_id

#Add ribosomal proteins that have been missed
fDat$Ribosomal <- fDat$Ribosomal | grepl("Rpl",fDat$symbol)

# Mitochondrial
fDat$Mitochondrial <- fDat$id %in% mitoGenes$V1

# TranscriptionFactor
fDat$TranscriptionFactor <- fDat$id %in% tfCheck$ensembl_gene_id

# Keep for hvg

fDat$KeepForHvg <- !(fDat$Ribosomal | fDat$Mitochondrial)


# Save data
stopifnot(identical(rownames(fDat),as.character(rownames(cDat))) & identical(colnames(cDat),pDat$barcode))
DataList <- list("phenoData"=pDat, "featureData"=fDat, "counts"=cDat)
saveRDS(DataList,file="../data/combined_Robjects/ExpressionList.rds")
