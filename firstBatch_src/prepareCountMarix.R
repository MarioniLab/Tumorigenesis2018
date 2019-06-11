
## This script extracts a cleaned count matrix from the 10X output
## The script is largely based on:
## https://github.com/MarioniLab/DropletUtils

# First remove reads derived from index swapping
library(DropletUtils)
samples <-  c("../../CellRangerAnalysis/SampleA11/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleB11/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleC11/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleD11/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleE11/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleF11/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleG11/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleH11/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleA12/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleB12/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleC12/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleD12/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleE12/outs/molecule_info.h5",
	      "../../CellRangerAnalysis/SampleF12/outs/molecule_info.h5")

bc.length <- 16
get.swapped <- FALSE
min.frac <- 0.9
cleanedCounts <- swappedDrops(samples=samples,
		     barcode.length=bc.length,
		     get.swapped=get.swapped,
		     min.frac=min.frac)
cleanedCounts <- cleanedCounts[[1]]
names(cleanedCounts) <- c("A11","B11","C11","D11","E11","F11","G11","H11","A12","B12","C12","D12","E12","F12")

finalCounts <- list()
# bpparam <- MulticoreParam(workers=3)

for (i in seq_along(cleanedCounts)) {
    m <- cleanedCounts[[i]]

    # test for empty droplets
    out <- emptyDrops(m)#,BPPARAM=bpparam)

    # set FDR threshold
    nonempty <- out$FDR < 0.01

    #remove NAs from NA FDR
    nonempty[is.na(nonempty)] <- FALSE

    # cant use logic vector to subset dgCMatrix?
    nonempty <- rownames(out[nonempty,])

    # Clean matrix
    m.clean <- m[,nonempty]
    finalCounts[[i]] <- m.clean
}

names(finalCounts) <- names(cleanedCounts)

m.out <- do.call(cbind,finalCounts)
ncells <- sapply(finalCounts,ncol)
smplnames <- rep(names(finalCounts),times=ncells)
colnames(m.out) <- paste(colnames(m.out),smplnames,sep="-")
saveRDS(m.out,"../data/firstBatch_Robjects/CountMatrix.rds")
