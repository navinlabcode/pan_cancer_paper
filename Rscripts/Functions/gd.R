#' Calculation and normalization of the geographic diversification (GD) metric
#' 
#' @param seg A segmented ratio matrix
#' @param times random sampling times
#' @param threads  number of CPU cores used for computation
#' 
#' @author GD method adapted from Wu et al. 2022, Cell Genomics
#' 

#Generate a distance matrix based on phylogenetic trees
sc_tree_dist<-function(df, dist_method='manhattan', tree_method='nj'){
  library(ape)
  #Calculate distance
  dist<-dist(df, method = dist_method)
  #Generating the tree using neighbor joining
  if (tree_method=='nj'){
    tree<-nj(dist)
  }
  #Or with Minimum evolution
  else{
    tree<-fastme.bal(dist)
  }
  return(cophenetic(tree))
}

#Based on the distance matrix, calculate GD
genomic_gd<-function(dist, k=10){
  #Order neighbors for each cell
  neighbors<-apply(dist, 2, order)
  neighbors<-sapply(1:ncol(neighbors), function(x){neighbors[,x][neighbors[,x]!=x]})
  #Spatial information
  sector<-str_extract(rownames(dist), "(S\\d+){1}")
  #Count number of cell pairs with K nearest neighbors that are in the same sector
  match<-0
  for(i in 1:ncol(neighbors)){
    match<-match+length(which(sector[neighbors[,i][1:k]]==sector[i]))
  }
  gd<-match/(k*ncol(neighbors))
  no_sector<-length(unique(sector))
  results<-list(gd=gd, no_sector = no_sector)
  return(results)
}

#Downsampling to normalize GD
gd_downsampling<-function(seg, times=100, threads=10){
  #Get sector information
  sector<-str_extract(rownames(seg), "(S\\d+){1}")
  #Random sample 4 sectors
  sectors_combn<-combn(unique(sector), 4)
  sectors_combn<-t(sectors_combn)
  set.seed(1245)
  sample<-sample.int(nrow(sectors_combn), size = times, replace = TRUE)
  #Calculate mean GD for all random samples
  results<-parallel::mclapply(sample, function(x) genomic_gd(dist=sc_tree_dist(df=seg[sector %in% sectors_combn[x,],])), 
                              mc.cores = threads)
  gd_mean<-mean(sapply(results, function(x) x[[1]]))
  #Make sure 4 sectors were sampled
  no_sector<-mean(sapply(results, function(x) x[[2]]))
  results<-list(gd=gd_mean, no_sector=no_sector)
  return(results)
}
