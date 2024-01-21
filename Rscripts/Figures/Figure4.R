#Get subclone colors
get_colors <- function() {
  hues1 <- paletteer::paletteer_d("ggsci::default_igv",
                                  n = 28) %>% unclass()
  
  
  colors <- setNames(hues1, paste0("c", 1:28))
  
  return(colors)
  
}

#Create gene annotations
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

#Chromosome annotations
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

#CN colors
wgd_cn_colors<-RColorBrewer::brewer.pal(11,  name = 'RdBu')[c(1:3,5,6,8,10)]
names(wgd_cn_colors)<-6:0

#Function to plot heatmaps with ploidy annotations
library(ComplexHeatmap)
plot_heatmap_int_wgd<-function(int, int.consensus, subclone_ploidy_long, ploidy_trunc=6, umap_data, colors=get_colors(), genes){
  subclone<-paste0('c', umap_data$cluster)
  subclone<-factor(subclone, levels = unique(subclone))
  #int.consesus is the consensus copy number profile for subclones
  rownames(int.consensus)<-paste0('c', rownames(int.consensus))
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
  #hclust to reorder all the cells
  hclust<-hclust(amap::Dist(int, method = 'manhattan', nbproc = 40), method = 'ward.D2')
  hclust_order<-hclust$labels[hclust$order]
  #int is the clustered copy number profile matrix
  int<-split(int, subclone)
  int<-lapply(int, function(x){x<-x[rownames(x)[order(match(rownames(x), hclust_order))],]})
  int<-int[tree_tips_order]
  names(int)<-NULL
  int<-do.call(rbind, int)
  #Truncate outlier copy numbers
  int[int>=ploidy_trunc]<-ploidy_trunc
  wgd<-rep('WGD', times=nrow(int))
  #Subclone ploidy long is the expanded subclonal consensus ploidy, just to help classify WGD and non-WGD subclones
  wgd[which(subclone_ploidy_long<2.5)]<-'nWGD'
  wgd<-wgd[order(match(rownames(umap_data), rownames(int)))]
  #ploidy here is the individual ploidy value for each single cell
  ploidy<-umap_data$ploidy
  #Truncate outlier ploidy values
  ploidy[ploidy>4.5]<-4.5
  ploidy[ploidy<1.5]<-1.5
  ploidy<-ploidy[order(match(rownames(umap_data), rownames(int)))]
  #Make gene annotations
  gene_annotation<-gene_anno_simple(gene_list = genes)
  ht<-Heatmap(as.matrix(int),top_annotation=chr_anno(),bottom_annotation = gene_annotation, name="Copy number", col=wgd_cn_colors, 
              row_title = paste0(nrow(int), " Single Cells"),row_title_gp = gpar(fontsize=20), row_title_rot  = 270,
              heatmap_legend_param = list(color_bar="discrete"), cluster_columns = FALSE,  show_row_names = FALSE,
              show_column_names = FALSE, use_raster = TRUE, raster_quality = 7, show_row_dend = F, cluster_rows = F, row_split = wgd)
  subclone<-subclone[order(match(subclone, tree_tips_order))]
  ha_cl<-HeatmapAnnotation(df = as.data.frame(subclone),col=list(subclone=colors) , which = "row", show_legend = T)
  myPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd")[-1])
  col_ploidy = circlize::colorRamp2(seq(1.5, 4.5, length=100), myPalette(100))
  ha_ploidy<-HeatmapAnnotation(df=as.data.frame(ploidy), col = list(ploidy=col_ploidy), which = 'row', show_legend = T)
  draw(ha_cl+ha_ploidy+ht, row_sub_title_side = "right" )
}

#Plot heatmap for G1
g1_genes<-c("DAXX", "GOPC", "HIF1A", "IDH1", "IDH2", "LZTR1", "MDM4", "PDGFRA", "PIK3CA",
            "PIK3R1", "ROS1", "SALL4", "STAG2", "TERT", "TP53") 
g1_int_clustered<-readRDS("Data/G1_int_clustered.rds")
g1_int_consensus<-readRDS("Data/G1_int_consensus.rds")
g1_umap<-readRDS("Data/G1_umap.rds")
g1_subclone_ploidy_long<-readRDS("Data/G1_subclone_ploidy_long.rds")
plot_heatmap_int_wgd(int=g1_int_clustered, int.consensus = g1_int_consensus, subclone_ploidy_long = g1_subclone_ploidy_long,
                     umap_data = g1_umap, genes = g1_genes)

#Plot heatmap for K7
k7_genes<-c("ARID1A", "BAP1", "CLTC", "HIF1A", "KDM5C", "KDM6A", "KMT2D", "MTOR", 
           "NF2", "PBRM1", "PTK6", "SETD2", "TFE3", "TSC1", "VHL", "TP53")  
k7_int_clustered<-readRDS("Data/K7_int_clustered.rds")
k7_int_consensus<-readRDS("Data/K7_int_consensus.rds")
k7_umap<-readRDS("Data/K7_umap.rds")
k7_subclone_ploidy_long<-readRDS("Data/K7_subclone_ploidy_long.rds")
plot_heatmap_int_wgd(int=k7_int_clustered, int.consensus = k7_int_consensus, subclone_ploidy_long = k7_subclone_ploidy_long,
                     umap_data = k7_umap, genes = k7_genes)

#Plot MEDICC2 trees
library(ggtree)
plot_medicc2_tree<-function(id, cluster){
  tree<-ape::read.tree(sprintf("Data/%s_medicc2_tree.new", id))
  #Count cell number to determine the tip size of the tree
  freq_df<-table(cluster)
  freq_df<-as.data.frame(freq_df)
  colnames(freq_df)<-c('taxa', 'cell_no')
  freq_df$taxa<-paste0('c', freq_df$taxa)
  freq_df$subclone<-freq_df$taxa
  ggtree::ggtree(tree) %<+% freq_df +  geom_tippoint(aes(size=cell_no, color=subclone))+scale_color_manual(values=get_colors())+
    scale_size(range = c(4, 10))
}

#G1
plot_medicc2_tree(id='G1', cluster = g1_umap$cluster)

#K7
plot_medicc2_tree(id='K7', cluster = g1_umap$cluster)

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

myPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd")[-1])
plot_umap_with_ploidy<-function(umap_data){
  #Re-order cells to leave outlier cells in the bottom of the umap plot
  umap_data<-split(umap_data, umap_data$cluster)
  umap_data<-lapply(umap_data, function(x) {x[order(x$ploidy, decreasing = ifelse(median(x$ploidy)<2.5, TRUE, FALSE)),]})
  umap_data<-do.call(rbind, umap_data)
  #Truncate outlier values
  umap_data$ploidy[umap_data$ploidy>4.5]<-4.5
  umap_data$ploidy[umap_data$ploidy<1.5]<-1.5
  my_theme <- list(ggplot2::theme(axis.title.x = element_text(colour = "gray28", size = 20), axis.text.x = element_blank(),
                                  axis.ticks.x = element_blank(),axis.title.y = element_text(colour = "gray28", size = 20), 
                                  axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(),
                                  legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 14),
                                  panel.border = element_rect(color = "black", fill = NA, size = 1)), 
                   xlab("UMAP1"), ylab("UMAP2"))
  print(ggplot(umap_data) + geom_point(aes(x = UMAP1, y = UMAP2, color = ploidy), size = 2)  + 
          scale_colour_gradientn(colors = myPalette(100), limits=c(1.5, 4.5)) + theme_classic() + my_theme)
}

#G1
plot_umap(g1_umap)
plot_umap_with_ploidy(g1_umap)

#K7
plot_umap(k7_umap)
plot_umap_with_ploidy(k7_umap)

#K8
k8_umap<-readRDS('Data/K8_umap.rds')
plot_umap(k8_umap)
plot_umap_with_ploidy(k8_umap)

#K9
k9_umap<-readRDS('Data/K9_umap.rds')
plot_umap(k9_umap)
plot_umap_with_ploidy(k9_umap)

#C1
c1_umap<-readRDS('Data/C1_umap.rds')
plot_umap(c1_umap)
plot_umap_with_ploidy(c1_umap)

#BR12
br12_umap<-readRDS('Data/BR12_umap.rds')
plot_umap(br12_umap)
plot_umap_with_ploidy(br12_umap)

#G2
g2_umap<-readRDS('Data/G2_umap.rds')
plot_umap(g2_umap)
plot_umap_with_ploidy(g2_umap)

#G3
g3_umap<-readRDS('Data/G3_umap.rds')
plot_umap(g3_umap)
plot_umap_with_ploidy(g3_umap)

#L7
l7_umap<-readRDS('Data/L7_umap.rds')
plot_umap(l7_umap)
plot_umap_with_ploidy(l7_umap)

#Plot summary
sWGD_summary<-readRDS("Data/Summary_sWGD.rds")
cancer_type<-c()
patient<-c()
cell_no<-c()
wgd_status<-c()
n_loss<-c()
for(i in 1:nrow(sWGD_summary)){
  cancer_type<-c(cancer_type, rep(sWGD_summary$cancer_type[i], 2))
  patient<-c(patient, rep(sWGD_summary$patient[i], 2))
  wgd_status<-c(wgd_status, 'nWGD', 'WGD')
  cell_no<-c(cell_no, sWGD_summary$cell_no_nWGD[i], sWGD_summary$cell_no_WGD[i])
  n_loss<-c(n_loss, sWGD_summary$n_loss_nWGD[i], sWGD_summary$n_loss_WGD[i])
}
wgd<-data.frame(matrix(ncol = 0, nrow = 36))
wgd<-cbind(cancer_type, patient, wgd_status, cell_no, n_loss, wgd)
wgd$wgd_status<-factor(wgd$wgd_status, levels = c('nWGD', 'WGD'))
summary<-readRDS("Data/summary.rds")
cancer_type_colors<-c('hotpink','slategray3', 'lawngreen', 'gold', 'thistle', 'mediumorchid2', 'dodgerblue2')
names(cancer_type_colors)<-unique(summary$cancer_type)

ggplot(wgd, aes(x = wgd_status, y = n_loss)) + 
  geom_boxplot(aes(fill = wgd_status), alpha=0.2, coef=0, show.legend = F, lwd=0.3) +
  geom_line(aes(group = patient), linewidth=0.2) + ylab('Number of chromosome losses') +
  geom_point(aes(color=cancer_type, size=cell_no))+ scale_size_continuous(range=c(1, 6))+
  scale_color_manual(values = cancer_type_colors)+ theme(axis.line = element_line(linetype = 1), axis.text.x=element_text(size=15), 
                                                         axis.title.y =element_text(size=15),axis.text.y=element_text(size=15),
                                                         panel.grid.minor = element_blank(), panel.border = element_blank(),panel.grid.major = element_blank(), 
                                                         panel.background = element_blank(), axis.title.x = element_blank())
#Statistical significance
wilcox.test(sWGD_summary$n_loss_nWGD, sWGD_summary$n_loss_WGD, paired = T, alternative = 'less')
