#Plot CNA classes
cna_class<-readRDS('Data/CNA_class.rds')
summary<-readRDS("Data/summary.rds")
cna_type<-rep(c('Clonal', 'Subclonal'), times=nrow(cna_class))
cancer_type<-c()
patient<-c()
count<-c()
for(i in 1:nrow(cna_class)){
  count<-c(count, cna_class$cCNA[i], cna_class$sCNA[i])
  patient<-c(patient, rep(as.character(rownames(cna_class)[i]), 2))
  cancer_type<-c(cancer_type, rep(summary$cancer_type[i], 2))
}

cna<-cbind(cancer_type, patient, cna_type, data.frame(count))

cna$patient<-factor(cna$patient, levels = rownames(cna_class)[order(cna_class$cCNA/(cna_class$cCNA+cna_class$sCNA))])
cna$cna_type<-factor(cna$cna_type, levels = c('Subclonal', 'Clonal'))

library(ggplot2)
ggplot(cna, aes(fill=cna_type, y=count, x=patient)) +
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c('deepskyblue', 'darkorange')) + facet_grid(~cancer_type, scales = 'free', space = 'free')+ 
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1, size = 7),  legend.title = element_blank(), 
        panel.spacing = unit(0.1, "lines"),  panel.grid.minor = element_blank(), panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.background = element_blank(),)+ xlab("")+ylab("Frequency")+labs(fill='CNA Type')+ 
  geom_hline(yintercept=0.05, linetype="dashed", color = "grey20") +scale_y_continuous(breaks=sort(c(seq(0, 1, length.out=5), 0.05)))

#Plot consensus heatmaps

#Get colors for subclones
get_colors <- function() {
  hues1 <- paletteer::paletteer_d("ggsci::default_igv",
                                  n = 28) %>% unclass()
  
  
  colors <- setNames(hues1, paste0("c", 1:28))
  
  return(colors)
  
}

#Gene annotation
gene_anno <- function(CNA_type, gene_list, 
                      bins_in_cna_pipeline=read.delim("Data/bins_in_cna_pipeline_bands.bed"), 
                      remove_Y = TRUE) {
  library(GenomicRanges)
  library(Homo.sapiens)
  suppressMessages(hg_genes_grange <- genes(Homo.sapiens, columns = "SYMBOL"))
  hg_genes_grange_df <- as.data.frame(hg_genes_grange) 
  hg_genes_grange_df$SYMBOL <- unlist(unname(as.vector(hg_genes_grange_df$SYMBOL)))
  gene_anno_df <- hg_genes_grange_df %>% dplyr::filter(SYMBOL %in% gene_list)
  names(gene_anno_df)[6] <- "gene"
  if (remove_Y == TRUE){
    bins_in_cna_pipeline <- bins_in_cna_pipeline %>% dplyr::filter(chr != "chrY")  
  }
  
  genes.annotation = character(length = nrow(bins_in_cna_pipeline))
  
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
  col_labels<-vector()
  for (i in 1:length(labels)){
    if (CNA_type$cna[which(genes.annotation==labels[i])]=='Clonal'){
      col_labels<-c(col_labels, 'darkorange')
    }
    else if (CNA_type$cna[which(genes.annotation==labels[i])]=='Subclonal'){
      col_labels<-c(col_labels, 'deepskyblue')
    }
    else{
      col_labels<-c(col_labels, "gray70")
    }
  }
  gene_anno_return <- HeatmapAnnotation(df=CNA_type, link = anno_mark(at = positions, labels = labels, labels_gp = gpar(fontsize=20, cex = 0.5, col =col_labels), side='bottom', link_width = unit(3, "mm")), 
                                        col = list(cna = c("Clonal" = "darkorange", "Subclonal" = "deepskyblue", "Neutral"="gray88")),
                                        which = 'column', show_annotation_name = F)
  return(gene_anno_return)
  
}

#Chromosome annotation
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
  # getting lengths for chr numbers annotation
  chr_rl_c <- c(1, cumsum(chr_lengths))
  # creating a data frame to calculate rowMeans
  chr_df <-  data.frame(a = chr_rl_c[1:length(chr_rl_c)-1],b= chr_rl_c[2:length(chr_rl_c)])
  chr_l_means <- round(rowMeans(chr_df))
  
  # creating the vector for chr number annotations
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

#Plot heatmap
plot_heatmap_seg_consensus<-function(seg.consensus, CNA_annotation, cluster, split, genes, colors=get_colors()){
  #Count cell number for each subclone
  cell_number<-janitor::tabyl(cluster) %>% dplyr::pull(n)
  names(cell_number)<-paste0('c', rownames(seg.consensus))
  #split is used to split two lineages
  names(split)<-paste0('c', rownames(seg.consensus))
  rownames(seg.consensus)<-paste0('c', rownames(seg.consensus))
  #Re-order cells
  #Minimum evolution trees to order subclones
  if(nrow(seg.consensus)>2){
    tree<-ape::ladderize(ape::fastme.bal(dist(seg.consensus, method = 'manhattan')))
    is_tip <- tree$edge[, 2] <= length(tree$tip.label)
    ordered_tips_index <- tree$edge[is_tip, 2]
    tree_tips_order <- tree$tip.label[ordered_tips_index] %>% rev()
  }
  else{
    tree_tips_order<-rownames(seg.consensus)
  }
  seg.consensus<-seg.consensus[tree_tips_order,]
  cell_number<-cell_number[tree_tips_order]
  split<-split[tree_tips_order]
  ht<-Heatmap(as.matrix(log2(seg.consensus)), top_annotation=chr_anno(), 
              bottom_annotation = gene_anno(CNA_type = data.frame(cna=CNA_annotation), gene_list = genes),
              col=circlize::colorRamp2(breaks = c(-1.5, 0, 1.5), c("dodgerblue3", "white", "firebrick3")), 
              name = 'Log2 (Ratio)', row_split = split, row_title = NULL,
              right_annotation= rowAnnotation(Number = anno_barplot(cell_number, gp = gpar(fill = colors[rownames(seg.consensus)]))),
              cluster_columns = FALSE,  show_row_names = FALSE, show_column_names = FALSE, use_raster = TRUE,
              raster_quality = 7, show_row_dend = F, cluster_rows = F)
  subclone<-factor(rownames(seg.consensus), levels = paste0('c', 1:nrow(seg.consensus)))
  ha_cl<-HeatmapAnnotation(df = as.data.frame(subclone), col=list(subclone=colors), which = "row", show_legend = T)
  myPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd")[-1])
  col_ploidy = circlize::colorRamp2(seq(1.5, 4.5, length=100), myPalette(100))
  draw(ha_cl+ht, row_sub_title_side = "right" )
}

#C2
colon_cancer_genes<-c('APC', 'TP53', 'KRAS', 'PIK3CA', 'MUC17', 'BRAF', 'FBXW7', 'SMAD4', 'RNF43', 'CHD4', 'ARID1A', 'SOX9')
c2_seg_consensus<-readRDS('Data/C2_seg_consensus.rds')
c2_umap<-readRDS('Data/C2_umap.rds')
#load pre-calculated CNA annotation 
c2_cna_annotation<-readRDS('Data/C2_CNA_annotation.rds')
plot_heatmap_seg_consensus(seg.consensus = c2_seg_consensus, cluster = c2_umap$cluster,
                           genes = colon_cancer_genes, CNA_annotation = c2_cna_annotation,
                           split = c(1, rep(2, 6)))

#C3
c3_seg_consensus<-readRDS('Data/C3_seg_consensus.rds')
c3_umap<-readRDS('Data/C3_umap.rds')
c3_cna_annotation<-readRDS('Data/C3_CNA_annotation.rds')
plot_heatmap_seg_consensus(seg.consensus = c3_seg_consensus, cluster = c3_umap$cluster,
                           genes = colon_cancer_genes, CNA_annotation = c3_cna_annotation,
                           split = c(rep(1, times=9), rep(2, times=9)))

#L8
l8_seg_consensus<-readRDS('Data/L8_seg_consensus.rds')
l8_umap<-readRDS('Data/L8_umap.rds')
l8_cna_annotation<-readRDS('Data/L8_CNA_annotation.rds')
l8_genes<-c("ALK", "CD74", "LRIG3", "SOX2", "TP53", "KRAS", "EGFR")
plot_heatmap_seg_consensus(seg.consensus = l8_seg_consensus, cluster = l8_umap$cluster,
                           genes = l8_genes, CNA_annotation = l8_cna_annotation,
                           split = c(rep(1, times=8), rep(2, times=3)))

#Plot number of mutations shared

#Load mutation data for L1 and L2 for C2
c2_L1_snv<-readRDS('Data/C2_L1_SNV.rds')
c2_L2_snv<-readRDS('Data/C2_L2_SNV.rds')

L1_mutations<-apply(c2_L1_snv, 1, function(x){paste0(x[1:5], collapse = '')})
L2_mutations<-apply(c2_L2_snv, 1, function(x){paste0(x[1:5], collapse = '')})

library(VennDiagram)
venn.plot1<-draw.pairwise.venn(area1 = length(L1_mutations), area2 = length(L2_mutations),
                                 cross.area = length(intersect(L1_mutations, L2_mutations)), 
                                 category = c("Lineage 1", "Lineage 2"), fill = c("orange", "blue"), rotation.degree = 180)

#C3
c3_L1_snv<-readRDS('Data/C3_L1_SNV.rds')
c3_L2_snv<-readRDS('Data/C3_L2_SNV.rds')

L1_mutations<-apply(c3_L1_snv, 1, function(x){paste0(x[1:5], collapse = '')})
L2_mutations<-apply(c3_L2_snv, 1, function(x){paste0(x[1:5], collapse = '')})

venn.plot2<-draw.pairwise.venn(area1 = length(L1_mutations), area2 = length(L2_mutations),
                               cross.area = length(intersect(L1_mutations, L2_mutations)), 
                               category = c("Lineage 1", "Lineage 2"), fill = c("orange", "blue"), rotation.degree = 180)

#L8
l8_L1_snv<-readRDS('Data/L8_L1_SNV.rds')
l8_L2_snv<-readRDS('Data/L8_L2_SNV.rds')

L1_mutations<-apply(l8_L1_snv, 1, function(x){paste0(x[1:5], collapse = '')})
L2_mutations<-apply(l8_L2_snv, 1, function(x){paste0(x[1:5], collapse = '')})

venn.plot3<-draw.pairwise.venn(area1 = length(L1_mutations), area2 = length(L2_mutations),
                               cross.area = length(intersect(L1_mutations, L2_mutations)), 
                               category = c("Lineage 1", "Lineage 2"), fill = c("orange", "blue"), rotation.degree = 180)
