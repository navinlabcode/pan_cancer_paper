#Read GD summary data
gd_summary<-readRDS('Data/GD_summary.rds')

#Plot GD by sector no
cancer_type_colors<-c('hotpink','slategray3', 'lawngreen', 'gold', 'thistle', 'mediumorchid2', 'dodgerblue2')
names(cancer_type_colors)<-unique(gd_summary$cancer_type)
ggplot(gd_summary[gd_summary$no_sector %in% c(4,6,8,12),], aes(y=gd, x=no_sector))+geom_boxplot(aes(x=no_sector), show.legend = FALSE, outlier.shape = NA)+theme_bw()+
  xlab("Number of sectors")+ylab("Normalized GD")+ theme(axis.line = element_line(linetype = 1), axis.text.x=element_text(size=15), 
                                                         axis.title.y =element_text(size=15), axis.title.x  =element_text(size=15),
                                                         axis.text.y=element_text(size=15), panel.grid.minor = element_blank(), 
                                                         panel.border = element_blank(),panel.grid.major = element_blank(),
                                                         panel.background = element_blank(), legend.position = 'none')+ 
  geom_point(position = position_jitter(seed=17, width = 0.2, height = 0.001), aes(color=cancer_type), size=2, alpha=1) +
  scale_color_manual(values = cancer_type_colors)

#Plot GD by cancer type
ggplot(gd_summary, aes(y= gd, x=cancer_type)) +geom_boxplot(aes(x=reorder(cancer_type, gd, FUN=median)),
                                                                      show.legend = FALSE, outlier.shape = NA)+theme_bw()+
  xlab("Cancer type")+ylab("Normalized GD")+ theme(axis.line = element_line(linetype = 1), axis.text.x=element_text(size=15), 
                                                   axis.title.y =element_text(size=15), axis.title.x  =element_text(size=15),
                                                   axis.text.y=element_text(size=15), panel.grid.minor = element_blank(), 
                                                   panel.border = element_blank(),panel.grid.major = element_blank(),
                                                   panel.background = element_blank() , legend.position = 'none')+ 
  scale_y_continuous(breaks = seq(0, 1, by=0.1))+
  geom_point(position = position_jitter(seed=17, width = 0.2, height = 0.001), aes(color=cancer_type), size=2, alpha=1)+
  scale_color_manual(values = cancer_type_colors)

#Plot GD with MPD
summary<-readRDS('Data/summary.rds')
gd_summary<-cbind(gd_summary, sapply(rownames(gd_summary), function(x) summary$MPD[summary$patient==x]))
colnames(gd_summary)[ncol(gd_summary)]<-'MPD'
ggplot(gd_summary)+geom_point(aes(x=MPD, y=gd, color=cancer_type))+
  geom_smooth(aes(x=MPD, y=gd),method='lm',color='black', formula= y~x)+
  ylab('Normalized GD')+xlab('MPD')+scale_color_manual(values=cancer_type_colors)+theme_classic()

#Plot heatmaps for two representative cases
get_colors <- function() {
  hues1 <- paletteer::paletteer_d("ggsci::default_igv",
                                  n = 28) %>% unclass()
  
  
  colors <- setNames(hues1, paste0("c", 1:28))
  
  return(colors)
  
}

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

gene_anno_simple <- function(gene_list,bins_in_cna_pipeline=read.delim("Data/bins_in_cna_pipeline_bands.bed"), remove_Y = TRUE) {
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
  gene_anno_return <- HeatmapAnnotation(link = anno_mark(at = positions, labels = labels, labels_gp = gpar(fontsize=20, cex = 0.5), side='bottom', link_width = unit(3, "mm")), col = list(CNA_annotation = c("Clonal" = "dimgrey", "Subclonal" = "gray88")),
                                        which = 'column', show_legend = F, show_annotation_name = F)
  return(gene_anno_return)
  
}

plot_heatmap_seg_consensus<-function(seg.consensus, cluster, genes, colors=get_colors()){
  #Count cell number for each subclone
  cell_number<-janitor::tabyl(cluster) %>% dplyr::pull(n)
  names(cell_number)<-paste0('c', rownames(seg.consensus))
  rownames(seg.consensus)<-paste0('c', rownames(seg.consensus))
  #Annotate genes
  gene_annotation<-gene_anno_simple(gene_list = genes)
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
  ht<-Heatmap(as.matrix(log2(seg.consensus)), top_annotation=chr_anno(), bottom_annotation = gene_annotation,
              col=circlize::colorRamp2(breaks = c(-1.5, 0, 1.5), c("dodgerblue3", "white", "firebrick3")), 
              name = 'Log2 (Ratio)', row_title = NULL,
              cluster_columns = FALSE,  show_row_names = FALSE, show_column_names = FALSE, use_raster = TRUE,
              raster_quality = 7, show_row_dend = F, cluster_rows = F)
  subclone<-factor(rownames(seg.consensus), levels = paste0('c', 1:nrow(seg.consensus)))
  ha_cl<-HeatmapAnnotation(df = as.data.frame(subclone), col=list(subclone=colors), which = "row", show_legend = T)
  ha_cell_no  = HeatmapAnnotation(Number=anno_barplot(x = cell_number, gp=gpar(fill=colors[tree_tips_order])),which = 'row')
  draw(ha_cl+ht+ha_cell_no, row_sub_title_side = "right")
}

#L9
l9_umap<-readRDS('Data/L9_umap')
l9_seg_consensus<-readRDS('Data/L9_seg_consensus.rds')
lung_cancer_genes<-c('TP53', 'EGFR', 'KRAS', 'STK11', 'KEAP1', 'MET', 'ERBB2', 'PIK3CA', 'BRAF', 'PTEN', 'RBM10', 'ATM', 'PTPRD',
                               'APC', 'ARID1A')
plot_heatmap_seg_consensus(seg.consensus = l9_seg_consensus, cluster = l9_umap$cluster, genes=lung_cancer_genes)

#K9
k9_umap<-readRDS('Data/K9_umap.rds')
k9_seg_consensus<-readRDS('Data/K9_seg_consensus.rds')
kidney_cancer_genes<-c("ARID1A", "BAP1", "CLTC", "HIF1A", "KDM5C", "KDM6A", "KMT2D", "MTOR", "NF2", "PBRM1",
                       "PTK6", "SETD2", "TFE3", "TSC1", "VHL", "TP53")
plot_heatmap_seg_consensus(seg.consensus = k9_seg_consensus, cluster = k9_umap$cluster, genes=kidney_cancer_genes)

#Make scatter pies
count_subclone<-function(df, no_subclone){
  freq<-lapply(1:no_subclone, function(x) length(which(df$cluster==x)))
  freq<-unlist(freq)
  freq<-freq/sum(freq)
}

make_scatter_pie<-function(umap, x_axis, y_axis, col_subclone=get_colors()){
  library("scatterpie")
  sector<-str_extract(rownames(umap), "(S\\d+){1}")
  sector<-factor(sector, levels = paste0('S', 1:length(unique(sector))))
  umap<-cbind(umap, sector)
  umap<-umap[order(umap$sector),]
  umap_list<-split(umap, umap$sector)
  no_sector<-length(umap_list)
  #Count frequencies of subclones in each sector
  subclone_freq<-lapply(umap_list, count_subclone, no_subclone=max(unique(umap$cluster)))
  subclone_freq<-as.data.frame(do.call(rbind, subclone_freq))
  subclone_names<-as.character(paste0('c', 1:ncol(subclone_freq)))
  region<-c(1:no_sector)
  subclone_freq<-cbind(x_axis,y_axis,region, subclone_freq)
  colnames(subclone_freq)<-c('x','y','region', subclone_names)
  radius<-as.numeric(lapply(umap_list,nrow))
  radius<-0.3*sqrt(radius*no_sector/nrow(umap))
  subclone_freq$radius<-radius
  ggplot() + geom_scatterpie(aes(x=x , y=y, group=region, r=radius), subclone_freq, cols = subclone_names,color=NA)+
    coord_equal()+ scale_fill_manual(values=col_subclone, name='Subclone')+
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_blank(), 
          axis.ticks = element_blank(), panel.background = element_blank()) # +  geom_scatterpie_legend(subclone_freq$radius, x=0, y=0, n=2, labeller=function(x) round(nrow(umap)*x/(0.3*no_sector))) 
  
}

#L9
make_scatter_pie(umap = l9_umap, x_axis = c(rep(-0.7,4), rep(0.7,4)), y_axis = c(4:1, seq(4.5, 1.5)))

#K9
make_scatter_pie(umap = k9_umap, x_axis = c(rep(0, times=3), rep(1, times=3), rep(0, times=3), rep(1, times=3)), 
                 y_axis = c(6:4, 6:4, 3:1, 3:1))
