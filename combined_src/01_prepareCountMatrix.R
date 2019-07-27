# First Batch
# First remove reads derived from index swapping
library(DropletUtils)
samples <-  c("../data/CellRangerOutput/firstBatch/SampleA11_molecule_info.h5",
	      "../data/CellRangerOutput/firstBatch/SampleB11_molecule_info.h5",
	      "../data/CellRangerOutput/firstBatch/SampleC11_molecule_info.h5",
	      "../data/CellRangerOutput/firstBatch/SampleD11_molecule_info.h5",
	      "../data/CellRangerOutput/firstBatch/SampleE11_molecule_info.h5",
	      "../data/CellRangerOutput/firstBatch/SampleF11_molecule_info.h5",
	      "../data/CellRangerOutput/firstBatch/SampleG11_molecule_info.h5",
	      "../data/CellRangerOutput/firstBatch/SampleH11_molecule_info.h5",
	      "../data/CellRangerOutput/firstBatch/SampleA12_molecule_info.h5",
	      "../data/CellRangerOutput/firstBatch/SampleB12_molecule_info.h5",
	      "../data/CellRangerOutput/firstBatch/SampleC12_molecule_info.h5",
	      "../data/CellRangerOutput/firstBatch/SampleD12_molecule_info.h5",
	      "../data/CellRangerOutput/firstBatch/SampleE12_molecule_info.h5",
	      "../data/CellRangerOutput/firstBatch/SampleF12_molecule_info.h5")

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

m.batch1 <- do.call(cbind,finalCounts)
ncells <- sapply(finalCounts,ncol)
smplnames <- rep(names(finalCounts),times=ncells)
colnames(m.batch1) <- paste(colnames(m.batch1),smplnames,sep="-")

# Second Batch
# First remove reads derived from index swapping
library(DropletUtils)
samples <-  c("../data/CellRangerOutput/secondBatch/SampleA2_molecule_info.h5",
	      "../data/CellRangerOutput/secondBatch/SampleB2_molecule_info.h5",
	      "../data/CellRangerOutput/secondBatch/SampleC2_molecule_info.h5",
	      "../data/CellRangerOutput/secondBatch/SampleD2_molecule_info.h5",
	      "../data/CellRangerOutput/secondBatch/SampleE2_molecule_info.h5",
	      "../data/CellRangerOutput/secondBatch/SampleF2_molecule_info.h5",
	      "../data/CellRangerOutput/secondBatch/SampleG2_molecule_info.h5",
	      "../data/CellRangerOutput/secondBatch/SampleH2_molecule_info.h5",
	      "../data/CellRangerOutput/secondBatch/SampleC3_molecule_info.h5")

bc.length <- 16
get.swapped <- FALSE
min.frac <- 0.9
cleanedCounts <- swappedDrops(samples=samples,
		     barcode.length=bc.length,
		     get.swapped=get.swapped,
		     min.frac=min.frac)
cleanedCounts <- cleanedCounts[[1]]
names(cleanedCounts) <- c("A2","B2","C2","D2","E2","F2","G2","H2","C3")


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

m.batch2 <- do.call(cbind,finalCounts)
ncells <- sapply(finalCounts,ncol)
smplnames <- rep(names(finalCounts),times=ncells)
colnames(m.batch2) <- paste(colnames(m.batch2),smplnames,sep="-")

stopifnot(!any(colnames(m.batch1) %in% colnames(m.batch2)))
stopifnot(identical(rownames(m.batch1),rownames(m.batch2)))
m.out <- cbind(m.batch1,m.batch2)
saveRDS(m.out,"../data/combined_Robjects/CountMatrix.rds")
