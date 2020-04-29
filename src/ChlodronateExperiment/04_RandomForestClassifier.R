library(scran)
library(scater)
library(randomForest)

sce <- readRDS("../../data/Integrated/Robjects/SCE_combined_final.rds")
set.seed(42)
rD <- data.frame(rowData(readRDS("../../data/Tumorigenesis/Robjects/SCE.rds")))
keep <- rownames(sce) %in% rownames(rD)[rD$KeepForHvg]
markers <- findMarkers(sce,sce$CellTypesFinal,block=sce$Batch,subset.row=keep)

genes <- lapply(markers, function(x) rownames(x)[1:20])
genes <- unique(as.character(unlist(genes)))

# Split into train and test set
labls <- factor(sce$CellTypesFinal)
set.seed(300)
sam <- sample(1:ncol(sce),floor(0.8*ncol(sce)))
train.labls <- labls[sam]
train.data <- logcounts(sce[genes,sam])
test.labls <- labls[-sam]
test.data <- logcounts(sce[genes,-sam])

# Construct random forest
library(randomForest)
Tree <- randomForest(x=as.matrix(t(train.data)),y=train.labls,
             xtest=as.matrix(t(test.data)),ytest=test.labls,keep.forest=TRUE)

saveRDS(Tree, "../../data/Integrated/Robjects/Tree.rds")
# Confusion matrix
# conf <- Tree$confusion
# conf[order(conf[,"class.error"]),]


# Variable importance
# fp <- Tree
# varImpPlot(fp,
#        sort=T,
#        n.var=55,
#        main="Top 50 - Variable Importance")


# ---- Classify new data ----
# pred <- predict(object=Tree, newdata=t(m.sub)) 
