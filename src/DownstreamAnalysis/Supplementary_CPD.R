library(ggplot2)
dot_plot = function(selected_rows = NULL,
                    selected_columns = NULL,
                    filename = 'plot.pdf',
                    width = 8,
                    height = 10,
                    means_path = './means.txt',
                    pvalues_path = './pvalues.txt',
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    output_extension = '.pdf',
		    sve=FALSE,
		    onlySig=TRUE,
		    getRows=FALSE
){

  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)

  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]

  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }

  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }

  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]

  if (onlySig) {
      oneSig <- rowSums(sel_pval<0.1)>=1 & rowSums(sel_means>0.5)>=1
      sel_pval <- sel_pval[oneSig,]
      sel_means <- sel_means[oneSig,]
      selected_rows <- selected_rows[oneSig]
      
  }

  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  #pr[pr==0] = 1
  #plot.data = cbind(plot.data,log2(pr))
  plot.data = cbind(plot.data,pr)
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

  p <-  ggplot(plot.data,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=mean)) +
  scale_color_gradientn('Mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=10, colour = "black"),
        axis.text.x = element_text(size= 10, angle = 90, hjust = 1),
        axis.text.y = element_text(size=8, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))

  if (sve) {
  if (output_extension == '.pdf') {
      ggsave(p,filename, width = width, height = height, device = cairo_pdf, limitsize=F)
  }
  else {
      ggsave(p,filename, width = width, height = height, limitsize=F)
  }
  } else {
  if (getRows) {
	  return(selected_rows)
      } else {
      return(p)
  } }
}


# Epithelial

dr <- "../../data/Downstream/CellPhoneDB_Epithelial_out/"
all_pval <- read.table(paste0(dr,'4dG/pvalues.txt'), header=T, stringsAsFactors = F, sep="\t", comment.char = '', check.names=F)
selected_columns <- c("LpAv|Hs","Hs|LpAv")
# Chosen interactions

rws1 <- c("IGF2_IGF2R", "IGF2_IGF1R", "ICAM1_AREG", "EGFR_AREG",
	  "TNFSF11_TNFRSF11A", "WNT4_FZD1","NOTCH1_WNT4",
	  "PTN_PTPRS")
rws2 <- c("CD44_FGFR2", "TGFB1_TGFbeta receptor1", "TGFB1_TGFbeta receptor2",
	 "LIFR_LIF","FGF1_FGFR2", "TIMP1_FGFR2")
rowlist <- list("HS_LP"=rws1,
		"LP_HS"=rws2)

for (nm in names(rowlist)) {
    selected_rows <- rowlist[[nm]]
    means_separator <- '\t'
    pvalues_separator <- '\t'

    # Prepare Bin1
    pvalues_path <- paste0(dr,"bin1/pvalues.txt")
    means_path <- paste0(dr,"bin1/means.txt")
    tum1_pval <- read.table(pvalues_path, header=TRUE, stringsAsFactors = FALSE, sep=means_separator, comment.char = '', check.names=FALSE)
    tum1_means <- read.table(means_path, header=TRUE, stringsAsFactors = FALSE, sep=pvalues_separator, comment.char = '', check.names=FALSE)

    intr_pairs <- tum1_pval$interacting_pair
    tum1_pval <- tum1_pval[match(selected_rows, intr_pairs), selected_columns]
    tum1_means <- tum1_means[match(selected_rows, intr_pairs), selected_columns]

    # Prepare Bin4
    pvalues_path <- paste0(dr,"bin4/pvalues.txt")
    means_path <- paste0(dr,"bin4/means.txt")
    tum4_pval <- read.table(pvalues_path, header=TRUE, stringsAsFactors = FALSE, sep=means_separator, comment.char = '', check.names=FALSE)
    tum4_means <- read.table(means_path, header=TRUE, stringsAsFactors = FALSE, sep=pvalues_separator, comment.char = '', check.names=FALSE)

    intr_pairs <- tum4_pval$interacting_pair
    tum4_pval <- tum4_pval[match(selected_rows, intr_pairs), selected_columns]
    tum4_means <- tum4_means[match(selected_rows, intr_pairs), selected_columns]

    # Prepare Gestation
    pvalues_path <- paste0(dr,"9dG/pvalues.txt")
    means_path <- paste0(dr,"9dG/means.txt")
    gest_pval <- read.table(pvalues_path, header=TRUE, stringsAsFactors = FALSE, sep=means_separator, comment.char = '', check.names=FALSE)
    gest_means <- read.table(means_path, header=TRUE, stringsAsFactors = FALSE, sep=pvalues_separator, comment.char = '', check.names=FALSE)

    intr_pairs <- gest_pval$interacting_pair
    gest9_pval <- gest_pval[match(selected_rows, intr_pairs), selected_columns]
    gest9_means <- gest_means[match(selected_rows, intr_pairs), selected_columns]

    # Prepare Gestation
    pvalues_path <- paste0(dr,"14dG/pvalues.txt")
    means_path <- paste0(dr,"14dG/means.txt")
    gest_pval <- read.table(pvalues_path, header=TRUE, stringsAsFactors = FALSE, sep=means_separator, comment.char = '', check.names=FALSE)
    gest_means <- read.table(means_path, header=TRUE, stringsAsFactors = FALSE, sep=pvalues_separator, comment.char = '', check.names=FALSE)

    intr_pairs <- gest_pval$interacting_pair
    gest_pval <- gest_pval[match(selected_rows, intr_pairs), selected_columns]
    gest_means <- gest_means[match(selected_rows, intr_pairs), selected_columns]

    # Combine
    colnames(gest9_pval) <- paste0("G9-",colnames(gest_pval))
    colnames(gest9_means) <- paste0("G9-",colnames(gest_means))
    colnames(gest_pval) <- paste0("G14-",colnames(gest_pval))
    colnames(gest_means) <- paste0("G14-",colnames(gest_means))
    colnames(tum4_pval) <- paste0("T4-",colnames(tum4_pval))
    colnames(tum4_means) <- paste0("T4-",colnames(tum4_means))
    colnames(tum1_pval) <- paste0("T1-",colnames(tum1_pval))
    colnames(tum1_means) <- paste0("T1-",colnames(tum1_means))

    sel_pval <- cbind(gest9_pval,gest_pval,tum1_pval,tum4_pval)
    sel_means <- cbind(gest9_means,gest_means,tum1_means,tum4_means)

    df_names <- expand.grid(selected_rows, colnames(sel_pval))
    pval <- unlist(sel_pval)
    pval[pval==0] <- 0.0009

    plot.data <- cbind(df_names,pval)
    pr <- unlist(as.data.frame(sel_means))
      #pr[pr==0] = 1
      #plot.data = cbind(plot.data,log2(pr))
    plot.data = cbind(plot.data,pr)
    colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
    my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)

    lvls <- c("G9-LpAv|Hs","G9-Hs|LpAv",
	      "G14-LpAv|Hs","G14-Hs|LpAv",
	      "T1-LpAv|Hs","T1-Hs|LpAv",
	      "T4-LpAv|Hs","T4-Hs|LpAv")
    plot.data$clusters <- factor(plot.data$clusters,
				     levels=lvls)

    p <-  ggplot(plot.data,aes(x=clusters,y=pair)) +
	geom_point(aes(size=-log10(pvalue),color=mean)) +
	scale_color_gradient2('Mean (Molecule 1, Molecule 2)', high="blue") +
	theme_bw() +
	theme(panel.grid.minor = element_blank(),
		panel.grid.major=element_line(colour="grey80",size=0.1,linetype="dashed"),
		axis.text=element_text(size=10, colour = "black"),
		axis.text.x = element_text(size= 10, angle = 90, hjust=.5),
		axis.text.y = element_text(size=8, colour = "black"),
		axis.title=element_blank(),
		panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
    ggsave(paste0("../../data/figures/Supplementary/Interactions_",nm,".pdf"),width=6.5,height=4)
}
