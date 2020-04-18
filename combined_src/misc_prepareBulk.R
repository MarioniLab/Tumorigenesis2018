# This script is just to format the DMBA data so I don't have to do it all the time
# This is the sample annotation
pDat <- read.csv("../../LinTracing/LinTracing/DMBATimeline/data/results/SLX-18286.HGCKCDRXX.s_1.contents.csv")
pDat$Barcode <- gsub("-","_",pDat$Barcode)
pDat$Sample <- paste0("S",pDat$Sample)
pDat$Timepoint <- substring(pDat$Sample,1,2)
pDat$Timepoint <- gsub("S","T",pDat$Timepoint)
pDat$Tumor <- pDat$Sample %in% c("S5BT","S5CT","S6DT")
pDat <- pDat[,c("Barcode","Sample","Timepoint","Tumor")]

# Read in Counts
counts <- read.csv("../../LinTracing/LinTracing/DMBATimeline/data/results/CountMatrix.csv",row.names=1)
counts <- counts[,pDat$Barcode]
rownames(pDat) <- colnames(counts) <- pDat$Sample

# Move stats to PD
pDat$Unmapped <- t(counts)[,1]
pDat$MultiMap <- t(counts)[,2]
pDat$NoFeat <- t(counts)[,3]
pDat$Ambiguous <- t(counts)[,4]
counts <- counts[-c(1:4),]
pDat$TotalReads <- colSums(counts)

# Add Gene Annotation
library(biomaRt)
ensembl  <- useMart("ensembl",dataset="mmusculus_gene_ensembl")#,host="useast.ensembl.org") 
gene.data <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), mart=ensembl)
gene.data <- gene.data[!is.na(gene.data[,1]),]
gene.data.extra <- data.frame("ensembl_gene_id"=setdiff(rownames(counts),gene.data$ensembl_gene_id),
			      "external_gene_name"=NA)
gene.data <- rbind(gene.data,gene.data.extra)
gene.data$uniq <- scater::uniquifyFeatureNames(gene.data[,1],gene.data[,2])
# Love this there are genes that are not in ensembl but are in the alignment
rownames(gene.data) <- gene.data$ensembl_gene_id
gene.data <- gene.data[rownames(counts),]

out <- list("counts"=counts,
	    "pDat"=pDat,
	    "fDat"=gene.data)
saveRDS(out,"../data/Robjects/DMBA_Bulk.rds")
