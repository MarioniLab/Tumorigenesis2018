# function to perform differential correlation test in two samples
# given a singlecellexperiment object

TwoSampleDCARS = function(sce,
                          assay = "logcounts",
                          PPI,
                          twoSample,
                          twoSampleRefLevel = NULL,
                          split = NULL,
                          cor_method = "spearman",
                          minNonZeroSamples = 10,
                          niter = 1000,
                          ncores = 4,
                          parallel = TRUE) {
  # sce is a SingleCellExperiment object
  # assay is the assay slot to extract the data
  # PPI is either a two column matrix with the gene pairs to test, 
  # or a single gene ID for which all other genes are tested with as well.
  # twosample is a character vector with the name of the factor for which sce should be split for the two groups
  # niter is the number of permutation iterations to perform
  # split is a character vector with the name of the factor for which the tests should be separated, e.g. different celltypes. If NULL then everything is tested.
  # minNonZeroSamples is minimum number of observations that should have nonzero values for
  # both the two genes for both of the two samples
  # ncores is the number of cores if running in parallel
  
  # sce = exprs_filt
  # sce$tomato_fac = ifelse(sce$tomato, "tompos", "tomneg")
  # assay = "logcounts"
  # PPI = "T"
  # twoSample = "tomato_fac"
  # twoSampleRefLevel = "tompos"
  # minNonZeroSamples = 10
  # split = "celltype.mapped"
  # niter = 1000
  # parallel = TRUE
  # ncores = 4
  
  if (class(PPI) == "character" & length(PPI) == 1) {
    singleGene = TRUE
    if (!PPI %in% rownames(sce)) {
      stop("PPI given as single gene but not found in rownames of sce")
    }
    PPI_mat_all = cbind(PPI, rownames(sce))
    rownames(PPI_mat_all) <- apply(PPI_mat_all,1,paste0, collapse = "_")
  } else {
    singleGene = FALSE
    PPI_mat_all = PPI
  }
  
  
  if (is.null(split)) {
    sce$split <- "all"
  } else {
    sce$split <- colData(sce)[,split]
  }
  
  sce$twoSample <- colData(sce)[,twoSample]
  if (length(unique(sce$twoSample)) > 2) {
    stop("More than two groups given!")
  }
  if (is.null(twoSampleRefLevel)) {
    twoSampleRefLevel = sort(unique(sce$twosample))[1]
    message(paste0("treating ", twoSampleRefLevel, " as the first sample"))
  }
  sce$isRef = sce$twoSample == twoSampleRefLevel
  
  splitgroups = as.character(sort(unique(sce$split)))
  message(paste0("testing ", length(splitgroups), " groups separately"))
  
  # this is the dataframe containing all results
  result_DF = NULL
  
  for (splitgroup in splitgroups) {
    message(paste0("testing ", splitgroup, " now..."))
    
    PPI_mat = PPI_mat_all
    
    # subset to the particular split group
    keep = sce$split == splitgroup & !is.na(sce$split)
    dat = assay(sce, assay)[,keep]
    isRef = sce$isRef[keep]
    
    if (min(table(isRef)) < minNonZeroSamples) {
      message("Not enough cells in either group for testing")
      next
    }
    
    if (singleGene) {
    # can do PPI filtering for singleGene, not supported for PPI given as matrix
      
    bindat_first = 1*(dat[PPI_mat[,2],isRef, drop = FALSE]>0)
    bindat_singleGene_first = 1*(dat[PPI,isRef, drop = FALSE]>0)
    num_nonzero_first = as.numeric(bindat_first %*% t(bindat_singleGene_first))
    names(num_nonzero_first) <- rownames(PPI_mat)
    
    bindat_second = 1*(dat[PPI_mat[,2],!isRef, drop = FALSE]>0)
    bindat_singleGene_second = 1*(dat[PPI,!isRef, drop = FALSE]>0)
    num_nonzero_second = as.numeric(bindat_second %*% t(bindat_singleGene_second))
    names(num_nonzero_second) <- rownames(PPI_mat)
    
    PPI_mat <- PPI_mat[num_nonzero_first >= minNonZeroSamples & num_nonzero_second >= minNonZeroSamples, , drop = FALSE]
    PPI_mat <- PPI_mat[!PPI_mat[,1] == PPI_mat[,2], , drop = FALSE]
    
    if (nrow(PPI_mat) == 0) {
      message("No gene pairs pass initial expression filter")
      next
    }
    
    message("now calculating correlations")
    
    calculateCorsDiff = function(dat, PPI_mat, isRef, testSingle = FALSE) {
      
      dat_first = dat[PPI_mat[,2],isRef]
      dat_second = dat[PPI_mat[,2],!isRef]
      
      singleGene_first = dat[PPI,isRef, drop = FALSE]
      singleGene_second = dat[PPI,!isRef, drop = FALSE]
      
      cors_first_obs = cor(t(as.matrix(dat_first)), t(as.matrix(singleGene_first)), method = cor_method)[,1]
      names(cors_first_obs) <- rownames(PPI_mat)
      cors_second_obs = cor(t(as.matrix(dat_second)), t(as.matrix(singleGene_second)), method = cor_method)[,1]
      names(cors_second_obs) <- rownames(PPI_mat)
      
      corsdiff_obs = cors_first_obs - cors_second_obs
      
      res = list(
        cors_first_obs = cors_first_obs,
        cors_second_obs = cors_second_obs,
        corsdiff_obs = corsdiff_obs
      )
      
      if (testSingle) {
        
        numeric_first = as.numeric(singleGene_first)
        numeric_second = as.numeric(singleGene_second)
        
        nonzero_cortest_first = apply(dat_first,1,function(x){
          cor.test(x, numeric_first, method = cor_method)$p.value
        })
        
        nonzero_cortest_second = apply(dat_second,1,function(x){
          cor.test(x, numeric_second, method = cor_method)$p.value
        })
        
        res[["nonzero_cortest_first"]] <- nonzero_cortest_first
        res[["nonzero_cortest_second"]] <- nonzero_cortest_second
      }
      
      return(res)
      
    }
    
    corsList = calculateCorsDiff(dat, PPI_mat, isRef, testSingle = TRUE)
  
    # dat_first = dat[PPI_mat[,2],isRef]
    # dat_second = dat[PPI_mat[,2],!isRef]
    # 
    # singleGene_first = dat[PPI,isRef, drop = FALSE]
    # singleGene_second = dat[PPI,!isRef, drop = FALSE]
    # 
    # cors_first_obs = cor(t(as.matrix(dat_first)), t(as.matrix(singleGene_first)), method = cor_method)[,1]
    # names(cors_first_obs) <- rownames(PPI_mat)
    # cors_second_obs = cor(t(as.matrix(dat_second)), t(as.matrix(singleGene_second)), method = cor_method)[,1]
    # names(cors_second_obs) <- rownames(PPI_mat)
    # 
    # corsdiff_obs = cors_first_obs - cors_second_obs
    
    result_df = data.frame(
      splitgroup = splitgroup,
      gene1 = PPI_mat[,1],
      gene2 = PPI_mat[,2]
    )
    
    result_df$cors_first_obs = corsList[["cors_first_obs"]]
    result_df$nonzero_cortest_first_pval = corsList[["nonzero_cortest_first"]]
    result_df$cors_second_obs = corsList[["cors_second_obs"]]
    result_df$nonzero_cortest_second_pval = corsList[["nonzero_cortest_second"]]
    result_df$corsdiff_obs = corsList[["corsdiff_obs"]]
    result_df$num_nonzero_first = num_nonzero_first[rownames(PPI_mat)]
    result_df$num_nonzero_second = num_nonzero_second[rownames(PPI_mat)]
    
    message(paste0("Now testing ", dim(PPI_mat)[1], " pairs of genes..."))

    # calculate permuted test statistics
    # want a list of length ncores with index of permutation with total niter
    # nperms = base::split(seq_len(niter), rep(seq_len(ncores), length.out = niter))
    
    permuted_corsdiff_list = mclapply(seq_len(niter), function(perm) {
      system("echo \"testing\" >> mclapply.txt")
      # then run  watch "wc -l mclapply.txt"
      res = calculateCorsDiff(dat, PPI_mat, sample(isRef))[["corsdiff_obs"]]
      return(res)
    }, mc.cores = ncores)
    
    permuted_corsdiff = do.call(cbind, permuted_corsdiff_list)
    
    perm_pval = rowMeans(abs(permuted_corsdiff) >= abs(corsList[["corsdiff_obs"]]))
    
    perm_sd = apply(permuted_corsdiff,1,sd)
    norm_pval = 2*pnorm(abs(corsList[["corsdiff_obs"]]), 
                       mean = 0, sd = perm_sd, lower.tail = FALSE)
    
    result_df$niter = niter
    result_df$perm_pval = perm_pval
    result_df$perm_sd = perm_sd
    result_df$norm_pval = norm_pval
    result_df$normFDR = p.adjust(norm_pval, method = "BH")
    
    result_DF <- rbind(result_DF, result_df)
    message(paste0("done testing for ", splitgroup, "!"))
    }
  }
  return(result_DF)  
}

# out = TwoSampleDCARS(sce, PPI = "T", twoSample = "tomato_fac", twoSampleRefLevel = "tompos",
#                      split = "celltype.mapped")


addFisherZTransformTest = function(sce,
                                   out,
                                   assayName,
                                   splitgroupName,
                                   tomatoName = "tomato_fac",
                                   gene1 = "T") {
  # sce is the single cell experiment object
  # out is the existing output dataframe
  # assayName is the assay to pull from
  # splitgroupName is the name of the column to do the testing with
  # tomatoName is the name of the column for the tomato
  # gene1 is the first gene to test against
  
  # assayName = "binarised"
  # splitgroupName = "celltype.mapped.masked"
  # tomatoName = "tomato_fac"
  
  dat = assay(sce, assayName)
  splitgroups = colData(sce)[,splitgroupName]
  tom = colData(sce)[,tomatoName]
  
  splitgroup_values = as.character(unique(splitgroups))
  splitgroup_values <- splitgroup_values[!is.na(splitgroup_values)]
  gene2_values = rownames(dat)
  
  fisher_matrix = matrix(NA, nrow = length(gene2_values),
                         ncol = length(splitgroup_values),
                         dimnames = list(gene2_values, splitgroup_values))
  
  # test per splitgroup
  for (splitgroup in splitgroup_values) {
    # splitgroup = "Anterior Primitive Streak"
    print(splitgroup)
    
    subdat = dat[,splitgroups == splitgroup & !is.na(splitgroups)]
    subtom = tom[splitgroups == splitgroup & !is.na(splitgroups)]
    
    # test per gene
    for (gene2 in gene2_values) {
      # gene2 = "Mrpl15"
      print(gene2)
      
      fisher_matrix[gene2, splitgroup] <- DCARS::fisherZtransformTest(subdat, gene1, gene2, subtom)
      
    }
  }
  
  # fdr per splitgroup
  fisher_matrix_fdr = apply(fisher_matrix,2,p.adjust, method = "BH")
  
  # get in same dim as out
  fisher_long = apply(out,1,function(x) fisher_matrix[as.character(x["gene2"]),as.character(x["splitgroup"])])
  
  fisher_fdr_long = apply(out,1,function(x) fisher_matrix_fdr[as.character(x["gene2"]),as.character(x["splitgroup"])])
  
  out$fisherZ_pval <- fisher_long
  out$fisherZ_fdr <- fisher_fdr_long
  
  return(out)
}
