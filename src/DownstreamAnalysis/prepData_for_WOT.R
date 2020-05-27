library(scran)
sce <- readRDS("../../data/Integrated/Robjects/SCE_combined_final.rds")
pca.tum <- readRDS("../../data/Tumorigenesis/Robjects/CorrectedPCA.rds")
ttime <- read.csv("../../data/Integrated/Robjects/TumorTime.csv",)[,-1]
ttime$ptimeBin[grepl("TM",ttime$Condition)] <- 5

library(dplyr)
smry <- group_by(ttime,ptimeBin) %>%
    summarize(time=round(mean(ptime),0))

ttime <- left_join(ttime,smry)

pD.t <- data.frame(colData(sce[,sce$Experiment=="Tumorigenesis"]))
pD.t <- dplyr::right_join(pD.t,ttime)
colnames(pD.t)[colnames(pD.t)=="barcode"] <- "id" # for WOT
pca.tum <- pca.tum[pD.t$id,]
pD.t$ptime <- round(pD.t$ptime,1)

write.csv(pca.tum,"../../data/Integrated/Robjects/Tumor_Pca_for_WOT.txt",quote=FALSE,row.names=TRUE)
write.csv(pD.t,"../../data/Integrated/Robjects/Tumor_pD_for_WOT.txt",quote=FALSE,row.names=FALSE)

# Pregnancy
pca.preg <- readRDS("../../data/Pregnancy/Robjects/CorrectedPCA.rds")
pD <- data.frame(colData(sce[,sce$Experiment=="Pregnancy"]))
colnames(pD)[colnames(pD)=="barcode"] <- "id" # for WOT
pca.preg <- pca.preg[pD$id,]
pD$time <- as.numeric(as.character(plyr::mapvalues(pD$Condition,
			   c("CTRL","4.5dG","9.5dG","14.5dG"),
			   c(0,4.5,9.5,14.5))))

write.csv(pca.preg,"../../data/Integrated/Robjects/Preg_Pca_for_WOT.txt",quote=FALSE,row.names=TRUE)
write.csv(pD,"../../data/Integrated/Robjects/Preg_pD_for_WOT.txt",quote=FALSE,row.names=FALSE)
