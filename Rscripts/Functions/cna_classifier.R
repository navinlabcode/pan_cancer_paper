#' Classifying CNA events
#' 
#' @param em A consensus integer copy number event matrix for subclones
#' @param ploidy DNA ploidy values of all subclones
#' @param sex Sex of the patient
#' @param min_bin Minimum number of bins for a segment to be included
#' 
#' @author Hanghui Ye
#' 

classify_cna<-function(em, ploidy, sex='F', min_bin=5){
  em<-as.data.frame(em)
  #save all bins for annotation purpose
  em_all<-em
  anno<-rep('Neutral', times=nrow(em_all))
  #Only keep segments with at least certain number of bins
  seg_kept<-which(as.numeric(em$n.probes)>=min_bin)
  em<-em[seg_kept,]
  #Remove chromosome annotations and number of bins info
  df<-em[,-(1:3)]
  df<-apply(df, 2, as.numeric)
  #df1 will be used to classify events based on df
  df1<-df
  cna<-c()
  # We want to know the maximum number of genome doubling in the tumor
  # If it is a non-WGD tumor
  if(all(ploidy<2.5)==TRUE){
    wgd<-0
  }
  # If it is a WGD tumor
  else{
    #At least one-time WGD
    wgd<-1
    #find subclonal WGD tumors, normalize populations with fewer WGD
    if(round(max(ploidy)/min(ploidy))==2){
      #Subclonal WGD with only one-time WGD
      if(min(ploidy)<2.5){
        df[,which(ploidy<2.5)]<-2*df[,which(ploidy<2.5)]
      }
      #If in a subclonal WGD tumor, there is WGD in the low ploidy population, then the maximum number should be 2 
      else{
        wgd<-2
        df[,which(ploidy<4)]<-2*df[,which(ploidy<4)]
      }
    }
    #Clonal WGD tumors, if the ploidy is higher > 6N, we also think it has 2 times WGD
    else{
      if (max(ploidy)>6){
        wgd<-2
      }
    }
    
  }
  #Generate breakpoint info with chromosome number and arms
  bp_df<-df
  for(i in 1:nrow(bp_df)){
    for(j in 1:ncol(bp_df)){
      bp_df[i,j]<-paste0(em[i, 1], em[i,2], ' ', bp_df[i,j])
    }
  }
  bp_df<-apply(bp_df, 2, function(x){
    rle_lengths<-rle(x)$lengths
    breakpoints<-c()
    for(i in 1:length(rle_lengths)){
      breakpoints<-c(breakpoints, rep(paste0(cumsum(rle_lengths)[i]-rle_lengths[i]+1, '-', cumsum(rle_lengths)[i]), 
                                      times=rle_lengths[i]))
    }
    return(breakpoints)
  })
  #Call non-WGD CNA events, ground state is 2*2^n, where x is the total times of WGD
  if(!sex %in% c('F', 'M')){
    stop('Sex info is wrong')
  }
  if(sex=='F'){
    x<-2
  }
  else{
    x<-1
  }
  chr_x<-which(em[,1]== 23)
  #Autosome
  for(i in 1:(min(chr_x)-1)){
    for(j in 1:ncol(df)){
      if(df[i,j]==2*2^wgd){
        df1[i,j]<-'Neu'
      }
      else if(df[i,j]>2*2^wgd){
        df1[i,j]<-'Amp'
      }
      else{
        df1[i,j]<-'Del'
      }
    }
  }
  #Chromosome X
  for(i in chr_x){
    for(j in 1:ncol(df)){
      if(df[i,j]==x*2^wgd){
        df1[i,j]<-'Neu'
      }
      else if(df[i,j]>x*2^wgd){
        df1[i,j]<-'Amp'
      }
      else{
        df1[i,j]<-'Del'
      }
    }
  }
  #Combine event classification with breakpoint info
  for(i in 1:nrow(df1)){
    for(j in 1:ncol(df1)){
      df1[i,j]<-paste0(df1[i,j], bp_df[i,j])
    }
  }
  #Annotate events
  k<-1
  for(i in seg_kept){
    if(length(unique(str_detect(df1[k,], 'Neu')))==1 && unique(str_detect(df1[k,], 'Neu'))==TRUE){
      anno[i]<-'Neutral'
    }
    else if(length(unique(df1[k,]))==1){
      anno[i]<-'Clonal'
    }
    else{
      anno[i]<-'Subclonal'
    }
    k<-k+1
  }
  anno<-rep(anno, times=as.numeric(em_all$n.probes))
  df1<-as.data.frame(df1)
  #List all CNA events
  total_CNA<-unlist(df1)
  total_CNA<-unique(total_CNA)
  for(i in 1:length(total_CNA)){
    #Neutral CNA
    if(str_detect(total_CNA[i], 'Neu')==TRUE){
      total_CNA[i]<-'nCNA'
    }
    #Clonal CNA
    else if(length(unique(apply(df1, 2, function(x){total_CNA[i] %in% x})))==1){
      total_CNA[i]<-'cCNA'
    }
    #Subclonal CNA
    else{
      total_CNA[i]<-'sCNA'
    }
  }
  cna_count<-c(length(which(total_CNA=='sCNA')), length(which(total_CNA=='cCNA')), length(which(total_CNA=='nCNA')))
  names(cna_count)<-c('sCNA', 'cCNA', 'nCNA')
  results<-list(anno=anno, cna_count=cna_count)
  return(results)
}