#Generate color scheme
get_colors <- function() {
  hues1 <- paletteer::paletteer_d("ggsci::default_igv",
                                  n = 28) %>% unclass()
  
  
  colors <- setNames(hues1, paste0("c", 1:28))
  
  return(colors)
  
}

#Generate gene annotations
gene_anno <- function(CNA_type, gene_list, 
                      bins_in_cna_pipeline=read.delim("Data/bins_in_cna_pipeline_bands.bed"), 
                      remove_Y = TRUE) {
  library(GenomicRanges)
  library(Homo.sapiens)
  #Get genome coordinates for the genes
  suppressMessages(hg_genes_grange <- genes(Homo.sapiens, columns = "SYMBOL"))
  hg_genes_grange_df <- as.data.frame(hg_genes_grange) 
  hg_genes_grange_df$SYMBOL <- unlist(unname(as.vector(hg_genes_grange_df$SYMBOL)))
  gene_anno_df <- hg_genes_grange_df %>% dplyr::filter(SYMBOL %in% gene_list)
  names(gene_anno_df)[6] <- "gene"
  if (remove_Y == TRUE){
    bins_in_cna_pipeline <- bins_in_cna_pipeline %>% dplyr::filter(chr != "chrY")  
  }
  genes.annotation = character(length = nrow(bins_in_cna_pipeline))
  #Find the bin for each gene
  for (i in 1:nrow(gene_anno_df)) {
    for (j in 1:nrow(bins_in_cna_pipeline)) {
      if ((findInterval(
        gene_anno_df$start[i],
        c(bins_in_cna_pipeline$start[j], bins_in_cna_pipeline$end[j])
      ) == 1) &&
      (bins_in_cna_pipeline$chr[j] == gene_anno_df$seqnames[i]))  {
        genes.annotation[j] <- paste0(genes.annotation[j], " ", as.character(gene_anno_df$gene[i]))
        break
      }
      else if((findInterval(
        gene_anno_df$end[i],
        c(bins_in_cna_pipeline$start[j], bins_in_cna_pipeline$end[j])
      ) == 1) &&
      (bins_in_cna_pipeline$chr[j] == gene_anno_df$seqnames[i])){
        genes.annotation[j] <- paste0(genes.annotation[j], " ", as.character(gene_anno_df$gene[i]))
        break
      }
    }
  }
  genes.annotation[genes.annotation == ""] <- NA
  labels <- genes.annotation[!is.na(genes.annotation)]
  positions <- which(!is.na(genes.annotation))
  padding = unit.c(unit(2, "cm"), unit(1, "cm"),
                   unit(c(2, 1), "cm"))
  #Annotate clonal or subclonal
  col_labels<-vector()
  for (i in 1:length(labels)){
    if (CNA_type$CNA_annotation[which(genes.annotation==labels[i])]=='Clonal'){
      col_labels<-c(col_labels, 'gray30')
    }
    else{
      col_labels<-c(col_labels, 'gray70')
    }
  }
  gene_anno_return <- HeatmapAnnotation(df=CNA_type, link = anno_mark(at = positions, labels = labels, labels_gp = gpar(fontsize=20, cex = 0.5, col =col_labels), side='bottom', link_width = unit(3, "mm")), col = list(CNA_annotation = c("Clonal" = "dimgrey", "Subclonal" = "gray88")),
                                        which = 'column', show_legend = F, show_annotation_name = F)
  return(gene_anno_return)
  
}

#Generate chromosome annotations
chr_anno<-function(remove_Y = TRUE){
  if(remove_Y == TRUE){
    chr_lengths <- c(953, 1038,  875,  838,  782,  724,  646,  635,  480,  561,  568,  575,  431,  391,  332,  310,  311,  336,  220, 266,  149 , 139,  607)
    chr_binary <-( (1:23) %% 2 )+1
    chrom.names <- c(1:22,"X")
  }
  else {
    chr_lengths<-c(953, 1038,  875,  838,  782,  724,  646,  635,  480,  561,  568,  575,  431,  391,  332,  310,  311,  336,  220, 266,  149 , 139,  607, 38)
    chr_binary<-(1:24)%%2 +1
    chrom.names <- c(1:22,"X", 'Y')
  }
  Chr <- data.frame(Chr = rep.int(x = chr_binary, times = chr_lengths))
  # getting lengths for chromosome number annotation
  chr_rl_c <- c(1, cumsum(chr_lengths))
  # creating a data frame to calculate rowMeans
  chr_df <-  data.frame(a = chr_rl_c[1:length(chr_rl_c)-1],b= chr_rl_c[2:length(chr_rl_c)])
  chr_l_means <- round(rowMeans(chr_df))
  # creating the vector for chromosome number annotations
  v<-vector(mode='character', length = cumsum(chr_lengths)[length(chr_lengths)])
  v[chr_l_means] <- chrom.names
  v[is.na(v)] <- ""
  chr_bar <- HeatmapAnnotation(show_annotation_name = FALSE, chr_text = anno_text(v, 
                                                                                  gp = gpar(fontsize = 8)),
                               df = Chr ,
                               show_legend = FALSE,
                               which = "column",annotation_name_side = 'left',
                               col = list(Chr = c("1" = "gray88", "2" = "black"))
  )
  return(chr_bar)
}

#Plot integer copy number heatmap, int: single cell integer copy number matrix, int.consensus: consensus integer copy number matrix
plot_heatmap_int<-function(int, int.consensus, ploidy_trunc, cluster, colors=get_colors(), genes){
  subclone<-paste0('c', cluster)
  rownames(int.consensus)<-paste0('c', rownames(int.consensus))
  subclone<-factor(subclone, levels = unique(subclone))
  #Create CNA annotation in the bottom
  CNA_annotation<-c()
  for(i in 1:ncol(int.consensus)){
    if(length(unique(int.consensus[,i]))==1){
      CNA_annotation<-c(CNA_annotation, 'Clonal')
    }
    else{
      CNA_annotation<-c(CNA_annotation, 'Subclonal')
    }
  }
  CNA_annotation<-as.data.frame(CNA_annotation)
  clonality<-gene_anno(CNA_type = CNA_annotation , gene_list = genes)
  #Re-order cells
  #Minimum evolution trees to order subclones
  if(nrow(int.consensus)>2){
    tree<-ape::ladderize(ape::fastme.bal(dist(int.consensus, method = 'manhattan')))
    is_tip <- tree$edge[, 2] <= length(tree$tip.label)
    ordered_tips_index <- tree$edge[is_tip, 2]
    tree_tips_order <- tree$tip.label[ordered_tips_index] %>% rev()
  }
  else{
    tree_tips_order<-rownames(int.consensus)
  }
  #hclust to order all single cells
  hclust<-hclust(amap::Dist(int, method = 'manhattan', nbproc = 30), method = 'ward.D2')
  hclust_order<-hclust$labels[hclust$order]
  int<-split(int, subclone)
  #Integrate orders of subclones and single cells
  int<-lapply(int, function(x){x<-x[rownames(x)[order(match(rownames(x), hclust_order))],]})
  int<-int[tree_tips_order]
  int<-do.call(rbind, int)
  #Truncate outlier copy number values
  int[int>=ploidy_trunc]<-ploidy_trunc
  color_heat <- structure(
    pals::ocean.balance(length(0:ploidy_trunc)),
    names = 0:(ploidy_trunc)  
  )
  ht<-Heatmap(as.matrix(int),top_annotation=chr_anno(),bottom_annotation = clonality, name="Copy number", col=color_heat, 
              row_title = paste0(nrow(int), " Single Cells"),row_title_gp = gpar(fontsize=20), row_title_rot  = 270,
              heatmap_legend_param = list(color_bar="discrete"), cluster_columns = FALSE,  show_row_names = FALSE,
              show_column_names = FALSE, use_raster = TRUE,raster_quality = 7, show_row_dend = F, cluster_rows = F)
  subclone<-subclone[order(match(subclone, tree_tips_order))]
  ha_cl<-HeatmapAnnotation(df = as.data.frame(subclone),col=list(subclone=colors) , which = "row", show_legend = T)
  draw(ha_cl+ht, row_sub_title_side = "right" )
}

#Plot UMAP
plot_umap<-function(umap_data, color_cluster=get_colors()){
  umap_data$cluster<-paste0('c', umap_data$cluster)
  umap_data$cluster<-factor(umap_data$cluster, levels = unique(umap_data$cluster))
  my_theme <- list(ggplot2::theme(axis.title.x = element_text(colour = "gray28", size = 20), axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(),axis.title.y = element_text(colour = "gray28", size = 20), 
                                  axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(),
                                  legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 14),
                                  panel.border = element_rect(color = "black", fill = NA, size = 1)), 
                   xlab("UMAP1"), ylab("UMAP2"))
  print(ggplot(umap_data) + geom_point(aes(x = UMAP1, y = UMAP2, color = cluster), alpha = 1, size = 2)  + scale_color_manual(values = color_cluster) + theme_classic() + my_theme)
}