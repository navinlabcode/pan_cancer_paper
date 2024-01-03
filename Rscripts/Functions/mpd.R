#' Calculate mean pairwise phylogenetic distance based on single cell segmented copy number ratio profiles
#'
#' @param df A segmented ratio matrix
#' @param dist_method Distance calculation method
#' @param tree_method  Method to generate the phylogenetic tree (Neighbor joining or Minimum evolution)
#' 
#' @author Junke Wang, Hanghui Ye
#' 

mpd_calculator<-function(df, dist_method='manhattan', tree_method='nj'){
  library(ape)
  #Calculated distances first
  dist<-dist(df, method = dist_method)
  #Generate the tree
  # neighbor joining
  if (tree_method=='nj'){
    tree<-nj(dist)
  }
  # Minimum evolution
  else{
    tree<-fastme.bal(dist)
  }
  tree_complete<-tree
  #Pruning trees
  n<-length(tree$tip.label)
  tree$edge.length[sapply(1:n,function(x,y)   which(y==x),y=tree$edge[,2])] <- 0
  #Calculate MPD
  mpd <- mean(cophenetic(tree))
  results<-list(tree=tree_complete, mpd=mpd)
  #Return tree and MPD
  return(results)
}
