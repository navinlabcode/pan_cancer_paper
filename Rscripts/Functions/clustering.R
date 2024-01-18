#' UMAP and HDBSCAN Clustering on segmented copy number ratios, plotting to visualize
#'
#' @param inuber A segmented ratio matrix
#' @param outprefix Output file name
#' @param plot Output plot name
#'
#' @author Darlan Conterno Minussi, Hanghui Ye
#' @return
#' @export
#' @examples
#' 
#' 

#Setup
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rwantshue))
suppressPackageStartupMessages(library(ComplexHeatmap))

parser <- ArgumentParser()

parser$add_argument("-i", "--inuber",
                    help = "A segmented ratio matrix")

parser$add_argument("-o", "--outprefix", type = "character",
                    help = "Output file name")

parser$add_argument("-p", "--plot", type = "character",
                    help = "Output plot name")


#Create color schemes
get_colors <- function(cluster_order_info,
                       seed = 1253) {
  
  color_space <- list(
    c(0, 360),	# 312hue range [0,360]
    c(0, 100),		# 36chroma range [0,100]
    c(0, 100))		# lightness range [0,100]
  scheme <- iwanthue(seed = seed, force_init = TRUE)
  
  hues1 <- scheme$hex(
    n = length(unique(cluster_order_info)),
    force_mode = FALSE,
    quality = 50,
    color_space = color_space)
  
  colors_vec1 <- setNames(hues1, levels(as.factor(cluster_order_info)))
  
  colors_list <- list(colors1 = colors_vec1)
  
  return(colors_list)
  
}

#Plot UMAP
plot_umap<-function(umap_data){
  color_cluster<-get_colors(umap_data$cluster)
  my_theme <- list(ggplot2::theme(axis.title.x = element_text(colour = "gray28", size = 20), axis.text.x = element_blank(), 
                                  axis.ticks.x = element_blank(),axis.title.y = element_text(colour = "gray28", size = 20), 
                                  axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(),
                                  legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 14)), 
                   xlab("UMAP1"), ylab("UMAP2"))
  print(ggplot(umap_data) + geom_point(aes(x = UMAP1, y = UMAP2, color = as.factor(cluster)), alpha = 1, size = 1)  
        + scale_color_manual(values = color_cluster$colors1) + theme_classic() + my_theme)
  return(color_cluster)
}

#Create chromosome annotations
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
  chr_bar <- HeatmapAnnotation(chr_text = anno_text(v, 
                                                    gp = gpar(fontsize = 8)),
                               df = Chr ,
                               show_legend = FALSE,
                               which = "column",annotation_name_side = 'left',
                               col = list(Chr = c("1" = "gray88", "2" = "black"))
  )
  return(chr_bar)
}

#Plot copy number heatmap for segmented ratios
plot_heatmap<-function(copy_number_data, umap_data, col_cluster){
  cluster<-umap_data$cluster
  #Annotation for clusters
  ha_cl<-HeatmapAnnotation(df = as.data.frame(cluster),col=list(cluster=col_cluster), which = "row", show_legend = T)
  #Create the heatmap
  ht<-Heatmap(as.matrix(log2(copy_number_data)),top_annotation=chr_anno(remove_Y = !(ncol(copy_number_data)==12205)), name="Log2 (Ratio)",
              row_title = paste0(nrow(umap_data), " Single Cells"), row_title_gp = gpar(fontsize=20), row_title_rot  = 270 , 
              col=circlize::colorRamp2(breaks = c(-1.5, 0, 1.5), c("dodgerblue3", "white", "firebrick3")),
              cluster_columns = FALSE, cluster_rows =FALSE, show_row_names = FALSE,show_column_names = FALSE, use_raster = TRUE, raster_quality = 5)
  draw(ha_cl+ht, row_sub_title_side = "right", padding = unit.c(unit(1.2, "cm"), unit(0.1, "cm"), unit(c(0.1, 0.1), "cm")))
} 

#Clustering
clustering<-function(inuber, outprefix, umap_dist='manhattan', k=13, seed_no=135, logr=FALSE, remove_Y=TRUE, plot){
  #read the segmented ratio file
  df<-read.delim(inuber)
  df<-dplyr::select(df, -chrom, -chrompos, -abspos)
  df<-as.data.frame(t(df))
  #remove chromosome Y
  if (remove_Y==TRUE){
    df<-df[1:12167]
  }
  set.seed(seed=seed_no)
  #We can choose Use log normalization or not
  if(logr==TRUE){
  dat_umap <- uwot::umap(log(df), metric = umap_dist,min_dist = 0,n_neighbors = 40, approx_pow = TRUE)
  }
  else{
  dat_umap <- uwot::umap(df, metric = umap_dist,min_dist = 0,n_neighbors = 40, approx_pow = TRUE)
  }
  umap_df <- as.data.frame(dat_umap)
  umap_df<-dplyr::rename(umap_df, "UMAP1"="V1", "UMAP2"="V2")
  rownames(umap_df) <- rownames(df)
  #HDBSCAN clustering
  hdb <- dbscan::hdbscan(umap_df, minPts = k)
  #If there is only one cluster
  if(length(unique(hdb$cluster))==1){
    hdb$cluster<-rep(1, times=length(hdb$cluster))
  }
  #Find the nearest clusters of the 'outliers' called by HDBSCAN
  dist_umap <- dist(umap_df[,c(1:2)]) %>% as.matrix() %>% as.data.frame()
  if(0 %in% hdb$cluster){
    outliers<-which(hdb$cluster==0)
    for (i in 1:length(outliers)){
      dist_umap_sliced<-dist_umap[-outliers, outliers[i]]
      nearest_cell<-rownames(dist_umap)[-outliers][which.min(dist_umap_sliced)]
      hdb$cluster[outliers[i]]<-hdb$cluster[which(rownames(dist_umap)==nearest_cell)]
    }
  }
  umap_df$cluster<-hdb$cluster
  umap_df$rowname<-rownames(umap_df)
  umap_df<-dplyr::arrange(umap_df,umap_df$cluster)
  rownames(umap_df)<-umap_df$rowname
  umap_df<-dplyr::select(umap_df, -rowname)
  #Plot UMAP and heatmap
  pdf(file= sprintf("%s", plot), height = 10, width = 8)
  message('Plotting figures')
  color_cluster<-plot_umap(umap_data = umap_df)
  df<-df[rownames(umap_df),]
  plot_heatmap(copy_number_data=df, umap_data=umap_df, col_cluster=color_cluster$colors1)
  dev.off()
  write.table(umap_df, sprintf("%s", outprefix))
  message('Clustering done')
}

message("Reading data.")

args <- parser$parse_args()

i <- args$inuber
o <- args$outprefix
p <- args$plot

if (is.null(o)) stop("No output filename provided, use flag -o")
if (is.null(p)) stop("No output figure name provided, use flag -p")

if (file.access(i) == -1) {
  stop(sprintf("File (%s) does not exist.", m))
} else {
  inuber <- i
}

clustering(inuber = inuber, outprefix = o, plot = p)

