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
mergeCluster <- function(x, clusters, min.DE=20, maxRep=10, removeGenes=NULL, merge=TRUE, ...)
{
    library(scran)
    counter <- c(1:maxRep)
    if (!is.null(removeGenes)) {
	x <- x[!(rownames(x) %in% removeGenes),]
    }
    for (i in counter) {
	clust.vals <- levels(clusters)
	marker.list <- findMarkers(x,clusters,subset.row=rowMeans(x)>0.01,...)

	out <- matrix(nrow=length(clust.vals), ncol=length(clust.vals))
	colnames(out) <- clust.vals
	rownames(out) <- clust.vals
	for (cl in names(marker.list)) {
	    sb <- marker.list[[cl]]
	    sb <- as.matrix(sb[,grepl("logFC",colnames(sb))])
	    if(ncol(sb)==1) {
		colnames(sb) <- paste0("logFC.",clust.vals[clust.vals!=cl])
	    }
	    sb <- colSums(abs(sb)>=1)
	    names(sb) <- substr(names(sb),7,nchar(names(sb)))
	    sb[cl] <- 0
	    sb <- sb[match(colnames(out),names(sb))]
	    out[,cl] <- sb
	}

	if (merge) {
	#Compute new groups
	hc <- hclust(as.dist(out))
	min.height <- min(hc$height)
	if (min.height <= min.DE) {
	    print(paste0(i,". joining of clusters with ",min.height," DE genes"))
	    cs <- cutree(hc,h=min.height)
	    newgrps <- sapply(c(1:max(cs)), function(i) paste(names(cs)[cs==i],collapse="."))
	    csnew <- mapvalues(cs,c(1:max(cs)), newgrps)
	    clusters <- mapvalues(clusters,names(csnew),csnew)
	} else {
	    break
	}
    } else {
	break}
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
