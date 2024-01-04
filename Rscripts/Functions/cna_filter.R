#' KNN filtering to remove outliers, followed by euploid/normal cells filtering
#'
#' @param inseg a segmented ratio matrix from CBS segmentation
#' @param inbin a bincount matrix from the VarBin pipeline
#' @param k_nearest number of nearest neighbors
#' @param resolution threshold of correlation values
#' @param n_threads number of threads
#' @param outprefix1 output bincount file name after outliers filtering
#' @param outprefix2 output bincount file name after euploid cells filtering
#' @param plot1 output heatmap name for outliers filtering
#' @param plot2 output heatmap name for euploid cells filtering
#'
#' @author Darlan Conterno Minussi, Hanghui Ye
#' @return
#' @export
#'
#' @examples
#' 


# SETUP

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ComplexHeatmap))

parser <- ArgumentParser()

parser$add_argument("-s", "--inseg",
                    help = "Segmented ratio matrix from CBS")

parser$add_argument("-b", "--inbin",
                    help = "A bincount matrix from the VarBin pipeline")

parser$add_argument("-f", "--plot1", type = "character",
                    help = "Output plot name -- filtering outliers")

parser$add_argument("-a", "--outprefix1", type = "character",
                    help = "Output bin file name -- containing all good cells")

parser$add_argument("-i", "--plot2", type = "character",
                    help = "Output plot name -- removing normal cells")

parser$add_argument("-t", "--outprefix2", type = "character",
                    help = "Output bin file name -- containing only tumor cells")

parser$add_argument("-n", "--n_threads", type = "integer",
                    help = "Number of threads. Default = 20",
                    default = 20)

parser$add_argument("-k", "--k_nearest", type = "integer",
                    help = "Number of Nearest neighbors",
                    default = 5)

parser$add_argument("-r", "--resolution", type = "double",
                    help = "Filtering Resolution",
                    default = 0.9)




cna_filter<-function (inseg, inbin, k_nearest = 5, resolution = 0.9, n_threads=20, outprefix1, 
                      outprefix2, plot1, plot2) 
{ 
  #Filter outliers
  message("Filtering outliers.")
  seg <- read.delim(inseg)
  bin <- read.delim(inbin)
  seg <- seg[,-(1:3)]
  # Remove Y before calculating correlation
  seg_s<-seg[1:12167,]
  # correction to avoid correlations calculations with standard deviation zero
  zero_sd_idx <- which(apply(seg_s, 2, sd) == 0)
  if (length(zero_sd_idx) >= 1) {
    seg_s[1, zero_sd_idx] <- seg_s[1, zero_sd_idx] + 1e-3
  }
  # Calculating correlation matrix
  dst <- cor(seg_s)
  dst_knn_df <- apply(as.matrix(dst), 1, function(x) {
    mean(sort(x, decreasing = T)[2:(k_nearest + 1)])
  }) %>% enframe(name = "sample", value = "cor")
  dst_knn_df <- dst_knn_df %>% mutate(filtered = case_when(cor >= resolution ~ "kept", cor < resolution ~ "removed"))
  sample<-colnames(seg)
  cell_filter<-tibble(sample)
  cell_filter$filtered <- dst_knn_df$filtered
  # Plot heatmap -- Filtering outliers
  message("Plotting heatmap -- Filtering outliers.")
  filter_anno <- rowAnnotation(filtered = dplyr::pull(dst_knn_df, filtered), col = list(filtered = c(kept = "green2", removed = "firebrick3")))
  #generating top annotation
  chr_lengths<-c(953, 1038,  875,  838,  782,  724,  646,  635,  480,  561,  568,  575,  431,  
                 391,  332,  310,  311,  336,  220, 266,  149 , 139,  607, 38)
  chr_binary<-(1:24)%%2 +1
  chrom.names <- c(1:22,"X", 'Y')
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
                                                    gp = gpar(fontsize = 8)), df = Chr, 
                               show_legend = FALSE, which = "column",annotation_name_side = 'left',
                               col = list(Chr = c("1" = "gray88", "2" = "black")))
  ht1 <- Heatmap(log2(t(seg)), cluster_rows = function(x) {
    fastcluster::hclust(amap::Dist(x, method = "manhattan", nbproc = n_threads), method = "ward.D2")
    },col=circlize::colorRamp2(breaks = c(-1.5, 0, 1.5), c("dodgerblue3", "white", "firebrick3")),
    row_title = paste0(ncol(seg), " Single Cells"), row_title_gp = gpar(fontsize=20), row_title_rot  = 90,
    cluster_columns = FALSE, use_raster = TRUE, raster_device = "CairoPNG",
    top_annotation = chr_bar, border = TRUE, show_row_names = FALSE, 
    show_column_names = FALSE, row_split = pull(dst_knn_df, filtered), 
    left_annotation = filter_anno, heatmap_legend_param = list(title = "log2(segratio)"))
  # plot the heatmap
  pdf(file=sprintf("%s", plot1), height = 10, width = 8)
  print(ht1)
  dev.off()
  cell_kept<-cell_filter$sample[which(cell_filter$filtered=="kept")]
  bin_filtered<-cbind(bin[,1:3],bin[cell_kept])
  readr::write_tsv(bin_filtered, sprintf("%s", outprefix1))
  message("Filtering outliers -- Done.")
  #Detect and filter normal cells
  message('Filtering normal cells')
  #remove XY
  seg<-seg[cell_kept]
  seg_s<-seg[1:11560,]
  # Calculate the coefficient of variation
  cv <- apply(seg_s, 2, function(x){sd(x)/mean(x)})
  # Add a simulated normal cells dataset to boost
  set.seed(17)
  cv_simul <- rnorm(1000,
                    mean = 0.015,
                    sd = 0.01)
  names(cv_simul) <- paste0("simul", 1:length(cv_simul))
  cv <- c(cv_simul, cv)
  fit <- mixtools::normalmixEM(cv)
  euploid_resolution <- fit$mu[1] + 5 * fit$sigma[1]
  cv <- cv[!grepl("simul", names(cv))]
  cv_df <- tibble::enframe(cv,
                           name = "sample",
                           value = "CV")
  cv_df_low_cv <- cv_df %>% mutate(is_normal = case_when(CV > euploid_resolution ~ 'tumor', CV <= euploid_resolution ~ 'normal'))
  message("Plotting heatmap -- Removing normal cells")
  filter_anno <- rowAnnotation(filtered = pull(cv_df_low_cv, is_normal), col = list(filtered = c(tumor = "green2", normal = "firebrick3")))
  ht2 <- Heatmap(log2(t(seg)), cluster_rows = function(x) {
    fastcluster::hclust(amap::Dist(x, method = "manhattan", nbproc = n_threads), method = "ward.D2")
  },col=circlize::colorRamp2(breaks = c(-1.5, 0, 1.5), c("dodgerblue3", "white", "firebrick3")),
  row_title = paste0(ncol(seg), " Single Cells"), row_title_gp = gpar(fontsize=20), row_title_rot  = 90, 
  cluster_columns = FALSE, use_raster = TRUE, 
  raster_device = "CairoPNG", top_annotation = chr_bar, border = TRUE, show_row_names = FALSE,
  show_column_names = FALSE, row_split = pull(cv_df_low_cv, is_normal), left_annotation = filter_anno, 
  heatmap_legend_param = list(title = "log2(segratio)"))
  #Plot heatmap
  pdf(file=sprintf("%s", plot2), height = 10, width = 8)
  print(ht2)
  dev.off()
  tumor_cell<-cv_df_low_cv$sample[which(cv_df_low_cv$is_normal=="tumor")]
  bin_filtered_t<-cbind(bin[,1:3],bin[tumor_cell])
  readr::write_tsv(bin_filtered_t, sprintf("%s", outprefix2))
  message("Removing normal cells -- Done.")
  }

# running filtering

message("Reading data.")

args <- parser$parse_args()

s <- args$inseg
b <- args$inbin
o1 <- args$outprefix1
o2 <- args$outprefix2
p1 <- args$plot1
p2 <- args$plot2
r <- args$resolution
n <- args$n_threads
k <- args$k_nearest


if (is.null(o1)) stop("No output filename provided, use flag -o")
if (is.null(p1)) stop("No output plot name provided, use flag -p")
if (is.null(o2)) stop("No output filename provided, use flag -o")
if (is.null(p2)) stop("No output plot name provided, use flag -p")


if (file.access(s) == -1) {
  stop(sprintf("File (%s) does not exist.", s))
} 
if (file.access(b) == -1){
  stop(sprintf("File (%s) does not exist.", b))
}
if (file.access(s) != -1 && file.access(b) !=-1){
  inseg <- s
  inbin <- b
}

cna_filter(inseg = inseg, inbin = inbin, resolution = r, outprefix1 = o1, outprefix2 = o2, 
           k_nearest = k, n_threads = n, plot1 = p1, plot2 = p2)


