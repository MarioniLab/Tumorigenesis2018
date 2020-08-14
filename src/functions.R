subSample <- function(dataList, cell.filter=NULL, cell.number=5000, group=NULL) {
    require(dplyr)
    pD <- dataList[["phenoData"]]
    if(!is.null(cell.filter)) pD <- pD[cell.filter,]
    if(!is.null(group)) pD <- group_by_(pD,group)
    cells <- sample_n(pD,size=cell.number, replace=FALSE) %>% .$barcode
    out <- list()
    out[["phenoData"]] <- pD[pD$barcode %in% cells,]
    out[["featureData"]] <- dataList[["featureData"]]
    out[["counts"]] <- dataList[["counts"]][,out[["phenoData"]][,"barcode"]]
    return(out)
}

getHighVar <- function(m, supress.plot=FALSE, get.var.out=FALSE, keepGenes=NULL, ...) {
    # wrapper function to get highVar genes
    # Highly variable genes 
    var.des <- trendVar(m,method="loess",min.mean=0.01,parametric=TRUE, ...)
    var.out <- decomposeVar(m,var.des)
    if (!is.null(keepGenes)) {
	stopifnot(length(keepGenes)==nrow(var.out))
	var.out <- var.out[keepGenes,]
    }
    o <- order(var.out$mean)
    hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >=0.5),]
    out <- rownames(hvg.out)
    if (!supress.plot) {
	plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
	    ylab="Variance of log-expression")
	lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
	points(hvg.out$mean, hvg.out$total, pch=16, col="red")
    }
    if (get.var.out) {
	return(var.out)
    } else {
	return(out)
    }
}
# this is adapted from scran with additional option to add a min value
isOut <- function(metric, nmads = 5, type = c("both", "lower", "higher"), 
                      log = FALSE, subset = NULL, batch = NULL, min_diff = NA,
		      min.value=NA) {
    if (log) {
        metric <- log10(metric)
    }
    if (any(is.na(metric))) { 
        warning("missing values ignored during outlier detection")
    }

    if (!is.null(batch)) {
        N <- length(metric)
        if (length(batch) != N) { 
            stop("length of 'batch' must equal length of 'metric'")
        }

        # Coercing non-NULL subset into a logical vector.
        if (!is.null(subset)) { 
            new.subset <- logical(N)
            names(new.subset) <- names(metric)
            new.subset[subset] <- TRUE
            subset <- new.subset
        }
   
        # Computing QC metrics for each batch. 
        by.batch <- split(seq_len(N), batch)
        collected <- logical(N)
        for (b in by.batch) {
            collected[b] <- Recall(metric[b], nmads = nmads, type = type,
                                   log = FALSE, subset = subset[b], 
                                   batch = NULL, min_diff = min_diff)
        }
        return(collected)
    }

    # Computing median/MAD (possibly based on subset of the data).
    if (!is.null(subset)) {
        submetric <- metric[subset]
        if (length(submetric) == 0L) {
            warning("no observations remaining after subsetting")
        }
    } else {
        submetric <- metric
    }
    cur.med <- median(submetric, na.rm = TRUE)
    cur.mad <- mad(submetric, center = cur.med, na.rm = TRUE)

    diff.val <- max(min_diff, nmads * cur.mad, na.rm = TRUE)
    upper.limit <- cur.med + diff.val
    lower.limit <- max(cur.med - diff.val, min.value, na.rm = TRUE)
    print(paste0("Lower limit is ", lower.limit))
    print(paste0("Upper limit is ", upper.limit))
    
    type <- match.arg(type)
    if (type == "lower") {
        upper.limit <- Inf
    } else if (type == "higher") {
        lower.limit <- -Inf
    }

    return(metric < lower.limit | upper.limit < metric)
}
# function to identify outliers after random forest prediction, based on Antonio's method
# Adopted from Nils
identify_unclassified<- function(pred.labels,true.labels,train.data,pred.data,selected.genes) {
    pred.labels <- as.character(pred.labels)
    for(i in levels(true.labels)){
      cur_cells_train <- train.data[selected.genes,true.labels == i]
      cur_cells_pred <- pred.data[selected.genes, pred.labels == i]

      # Calculate the median pairwise dissimilarity between each pred cell and all tainings cells
      medians <- apply(cur_cells_pred, 2, function(n){
	cur_distances <- apply(cur_cells_train, 2, function(x){
	  sqrt((1-cor(n, x))/2)
	})
	median(cur_distances)
      })
      
      # Find maximum pairwise distance within the cluster
      dist.all <- as.dist(sqrt((1 - cor(cur_cells_train))/2))
      
      # Assign outlying cells as unclassified
      short_labels <- pred.labels[pred.labels == i]
      short_labels[as.numeric(medians) > max(dist.all)] <- "Unclassified"
      pred.labels[pred.labels == i] <- short_labels
    }
    return(as.factor(pred.labels))
}

#Cluster stuff
mergeCluster <- function(x, clusters, min.DE=10, maxRep=20, removeGenes=NULL, merge=TRUE, ...)
{
    # This function iteratively merges clusters with a number of DE genes less or equal than min.DE
    # Clusters is a factor corresponding with a length equal to ncol(x)
    # x is the log-transformed normalized gene expression matrix (genes x cells)
    # DE is defined as lfc>1 and log FDR <= -2 
    require(scran)
    counter <- c(1:maxRep)
    # if uninteresting features are to be removed
    if (!is.null(removeGenes)) {
	x <- x[!(rownames(x) %in% removeGenes),]
    }
    # Iterate through computing DE genes and merging clusters (if merge=TRUE)
    for (i in counter) {
	# Compute DE Genes
	clust.vals <- levels(clusters)
	marker.list <- findMarkers(x,clusters,subset.row=rowMeans(x)>0.01, lfc=1, full.stats=TRUE, ...)

	## Create Cluster x Cluster matrix with each entry corresponding to number of DE genes (log.FDR <=-2)
	# Setup matrix
	out <- matrix(nrow=length(clust.vals), ncol=length(clust.vals))
	colnames(out) <- clust.vals
	rownames(out) <- clust.vals
	# loop through each cluster and get DE genes against every other cluster
	for (cl in names(marker.list)) {
	    sb <- marker.list[[cl]]
	    sb <- data.frame(sb)
	    sb <- as.data.frame(sb[,grepl(".log.FDR",colnames(sb))]) # as.data.frame in case it's a single column
	    if(ncol(sb)==1) {
		colnames(sb) <- setdiff(names(marker.list),cl)
	    } else {
	    colnames(sb) <- gsub("stats.|.log.FDR","",colnames(sb))
	    }
	    sb <- colSums(as.matrix(sb) <= -2) # -2 is the log FDR threshold for DE
	    sb[cl] <- 0
	    sb <- sb[match(colnames(out),names(sb))]
	    out[,cl] <- sb
	}

	if (merge) {

	# If a test does not have enough d.f. then findMarkers will output NA
	# This usually happens if you use say a blocking factor and 
	# some clusters only appear in some level of that factor
	# In that case I don't have information to merge them
	if (any(is.na(out))) { 
	    print("A test did not have enough D.F.") 
	    out[is.na(out)] <- min.DE+1 
	}
	#Compute new groups
	if (any(is.na(out))) { print("A test did not have enough D.F.") }
	out[is.na(out)] <- min.DE+1 # This is to prevent hclust from crashing, NAs are produced when their are not enough d.f. for testing. In this case I do not wnat the clusters to be merged because I formally cant test for DE.
	hc <- hclust(as.dist(out))
	#get min number of DE genes
	min.height <- min(hc$height)
	if (min.height <= min.DE) {
	    print(paste0(i,". joining of clusters with ",min.height," DE genes"))
	    cs <- cutree(hc,h=min.height)
	    newgrps <- sapply(c(1:max(cs)), function(i) paste(names(cs)[cs==i],collapse="."))
	    csnew <- plyr::mapvalues(cs,c(1:max(cs)), newgrps)
	    clusters <- plyr::mapvalues(clusters,names(csnew),csnew)
	} else {
	    break
	}
    } else {
	break}

	if (length(unique(clusters))==1) {
	    break
	}
    }
	output <- list("Marker"=marker.list,
		       "NumberDE"=out,
		       "NewCluster"=clusters)
	return(output)
}

subCluster <- function(m, clusters, removeGenes=NULL, method="dynamic", ...) {
    require(cluster)
    require(dynamicTreeCut)
    # Define clustering function
    if (method!="dynamic") {
    func <- function(x,k) {
	dis <- as.dist(sqrt((1-cor(t(x)))/2))
	tree <- hclust(dis, method="ward.D2")
	cluster <- cutree(tree,k=k)
	out <- list()
	out$cluster <- cluster
	return(out)
    }
    } else {
	func <- function(x) {
	    dis <- as.dist(sqrt((1-cor(x))/2))
	    tree <- hclust(dis, method="ward.D2")
	    maxCoreScatter <- 0.3
	    minGap <- (1 - maxCoreScatter) * 3/4
	    cs <- cutreeDynamic(tree,distM=as.matrix(dis),deepSplit=0, minClusterSize=max(50,floor(ncol(x)/100)), maxCoreScatter=maxCoreScatter, minGap=minGap)
	    out <- list()
	    out$cluster <- cs
	    return(out)
	}
    }

    # Perform subclustering
    out <- NULL
    for (cl in levels(clusters)) {
	cells <- colnames(m)[clusters==cl]
	m.sub <- m[,cells]
	m.sub <- m.sub[rowMeans(m.sub) >0.01,]

	# Highly variable genes
	hvg <- getHighVar(m.sub,get.var.out=TRUE,supress.plot=TRUE)
	if (!is.null(removeGenes)) {
	hvg <- hvg[!(rownames(hvg) %in% removeGenes),]	
	}
	hvg <- rownames(hvg[order(hvg$bio, decreasing=TRUE),])[1:(nrow(hvg)/10)]

	x <- m.sub[hvg,]
	if (method!="dynamic") {
	gpas <- clusGap(t(x), func, K.max=3, B=2)
	k.opt <- maxSE(gpas$Tab[,"gap"],gpas$Tab[,"SE.sim"])
	subclust <- func(t(x),k=k.opt)
	} else {
	    subclust <- func(x)
	}
	tmp <- paste(cl,subclust$cluster,sep=".")
	#         tmp <- paste(cl,subclust,sep="-")
	names(tmp) <- cells
	out <- c(out,tmp)
    }
    out <- out[colnames(m)]
    return(as.factor(out))
}
bubblePlot <- function(m, markers, grps, cluster_col=TRUE, cluster_row=TRUE, angled=TRUE, zscores=FALSE) { 
    out <- data.frame(numeric(length(markers)))
    colnames(out) <- levels(factor(grps))[1]
    out.freq <- out
    for (cond in levels(factor(grps))) {
	expr <- rowMeans(m[markers,grps==cond])
	expr.freq <- rowMeans(m[markers,grps==cond]>0)
	colname <- cond
	out[,colname] <- expr
	out.freq[,colname] <- expr.freq
    }
    out <- as.matrix(out)
    out.freq <- as.matrix(out.freq)
    rownames(out.freq) <- rownames(out) <- markers
    if(zscores) {
	out <- t(scale(t(out))) 
    } else {
	out <- out/rowMax(out)
    }
    out.long <- melt(out,value.name="Mean")
    out.freq <- melt(out.freq,value.name="Frequency")
    out.long$Frequency <- out.freq$Frequency * 100

    if (cluster_row) {
    dis <- dist(out)
    hclst <- hclust(dis,method="ward.D2")
    lvlsVar1 <- rownames(out)[hclst$order]
    out.long$Var1 <- factor(out.long$Var1, levels=lvlsVar1)
    }
    if (cluster_col) {
    dis <- dist(t(out))
    hclst <- hclust(dis,method="ward.D2")
    lvlsVar2 <- colnames(out)[hclst$order]
    out.long$Var2 <- factor(out.long$Var2, levels=lvlsVar2)
    }
    p <- ggplot(out.long, aes(x=Var2, y=Var1, color=Mean, size=Frequency)) +
    geom_point() +
    scale_color_gradient2(low="grey90",high="black") + ##,mid="orange") +
   # scale_color_distiller(palette="Spectral") +
    scale_size(range=c(0,3)) +
    theme(panel.grid.major=element_line(colour="grey80",size=0.1,linetype="dashed"),
	  axis.text=element_text(size=7),
	  legend.position="bottom",
	  legend.direction="horizontal") +
    xlab("") +
    ylab("")
    if (angled) {
    p <- p + theme(axis.text.x=element_text(size=9, angle=45, hjust=1, vjust=1))
    }
    return(p)
}

quickGSE <- function(gens, pvals, ont="BP", method=c("GO", "Broad")) {
    # The sign selects the genes, all positive pvals between 0 and 0.01 will be regarded as significant
    require(org.Mm.eg.db)
    require(topGO)
    topDiffGenes <- function(allScore) { # topGO requirement..this is so useless..
	return(allScore < 0.01 & allScore >= 0)
    }
    alG <- pvals
    names(alG) <- gens


    if (grepl("GO",method)) {
	# prepare Data for topGO
	GO.data <- new("topGOdata", description="Lib GO",ontology=ont, allGenes=alG, 
		       annot=annFUN.org, mapping="org.Mm.eg.db",geneSelectionFun=topDiffGenes,
		       nodeSize=5, ID="symbol")
	# Fisher
	result.classic <- runTest(GO.data, statistic="fisher",algorithm="weight01")
	# KS I am just leaving this here in case I want to do this again
	#         test.stat <- new("classicScore", testStatistic=GOKSTest, name="KS tests")
	#         result.ks <- getSigGroups(GO.data,test.stat)
	# Output
	output <- GenTable(GO.data, Fisher=result.classic, #KS=result.ks, orderBy="KS",
			   orderBy="Fisher",topNodes=50, numChar=300)
	output$Term <- factor(output$Term, levels=unique(rev(output$Term)))
	return(output)
    } 

    if (method=="Broad") {
	require(reshape2)
	require(clusterProfiler)
	load("~/Cambridge/Tumorigenesis2018/data/misc/mouse_H_v5p2.rdata")
#	load("~/Cambridge/Tumorigenesis2018/data/misc/mouse_c2_v5p2.rdata")
#	load(url("http://bioinf.wehi.edu.au/software/MSigDB/mouse_c7_v5p2.rdata"))

	geneset <- melt(Mm.H)
#	geneset2 <- melt(Mm.c2)
#	geneset3 <- melt(Mm.c7)
#	geneset <- rbind(geneset1,geneset2)#,geneset3)

	colnames(geneset) <- c("ENTREZID","Geneset")
	ids <- unique(as.character(geneset$ENTREZID))
	res <- select(org.Mm.eg.db, keys=ids, column=c("SYMBOL"))
	geneset <- full_join(geneset,res)
	geneset <- geneset[!duplicated(geneset),]
	geneset <- geneset[,c(2:3)]
	colnames(geneset) <- c("TERM","GENE")
	out <- enricher(names(alG)[topDiffGenes(alG)],universe=names(alG),TERM2GENE=geneset,qvalueCutoff=0.01,minGSSize=10)
	out <- data.frame(out)
	out$ID <- factor(out$ID, levels=rev(unique(out$ID)))
	return(out)
    }
}

sumCountsAcrossCells  <- function(m,grps) {
    # Function to create pseudo-bulk samples from groups (e.g. clusters, samples) of cells
    grps <- factor(grps) #drop unused levels
    out <- data.frame(numeric(nrow(m)))
    colnames(out) <- levels(grps)[1]
    for (clust in levels(grps)) {
	expr <- Matrix::rowSums(m[,grps==clust])
	colname <- clust
	out[,colname] <- expr
    }
    rownames(out) <- rownames(m)
    return(out)
}

get_umap_graph <- function(umap_output) {
    require(Matrix)
    require(igraph)
    # extracts the graph information and outputs a weighted igraph
    dists.knn <- umap_output$knn[["distances"]]
    indx.knn <- umap_output$knn[["indexes"]]
    m.adj <- Matrix(0, nrow=nrow(indx.knn), ncol=nrow(indx.knn), sparse=TRUE) 
    rownames(m.adj) <- colnames(m.adj) <- rownames(indx.knn)

    for (i in seq_len(nrow(m.adj))) {
	m.adj[i,rownames(indx.knn)[indx.knn[i,]]] <- dists.knn[i,] 
    }

    igr <- graph_from_adjacency_matrix(m.adj, weighted=TRUE)
    return(igr)
}

theme_pub <- function(base_size=16) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
  + theme(plot.title = element_text(face = "bold"),
          text = element_text(size=18),
          panel.background = element_rect(colour = NA),
          plot.background = element_rect(colour = NA),
          panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1.2)),
          axis.title.y = element_text(angle=90,vjust =2),
          axis.title.x = element_text(vjust = -0.2),
          axis.text = element_text(size=rel(1)), 
          axis.line = element_line(colour="black",size=rel(1.5)),
          axis.ticks = element_line(),
          panel.grid.major = element_line(colour="#f0f0f0"),
          panel.grid.minor = element_blank(),
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.key.width= unit(1.0, "cm"),
          legend.key.height = unit(0.5, "cm"),
          legend.title = element_text(face="italic"),
          strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold")
  ))
}
