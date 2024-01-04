#' Merge levels to smooth population segmented ratios from multipcf
#'
#' @param inbin a bincount matrix from the VarBin pipeline
#' @param output output file name
#' @param popseg population segmentation file from multipcf
#' @param cores number of threads
#'
#' @return
#' @export
#'
#' @examples
#' 



# Setup
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(parallel))

parser <- ArgumentParser()

parser$add_argument("-b", "--bincounts",
                    help = "A bin count matrix resulting from the VarBin pipeline")

parser$add_argument("-o", "--output", type = "character",
                    help = "Output file name")

parser$add_argument("-p", "--popseg",
                    help = "A population segmentation file from multipcf.")

parser$add_argument("-c", "--cores", type = "integer",
                    help = "Number of threads. Default = 20",
                    default = 20)


# Make popseg from short to long file

# this function generates long files (one column per segment) using the short files
create_popseg_long <- function(popseg_file) {
  popseg_times <- popseg_file$n.probes
  popseg_clean <- dplyr::select(popseg_file, -c("chrom", "arm", "start.pos", "end.pos", "n.probes"))
  popseg_t <- as.data.frame(t(popseg_clean))
  popseg_long <- as.data.frame(apply(popseg_t, 1, function(m) { rep.int(m, popseg_times) }))
  popseg_long_t <- as.data.frame(t(popseg_long))  
  return(popseg_long_t)
}

# Merge Levels

#rkhorse function:
#Combining two levels with given medians in a given samples
#vecObs: vector of observed log2ratios
#vecPredNow: current vector of predicted values (median of the corresponding
#            segment
#mnNow:  current vector of the segment values (medians)
#mn1: median of the first segment to test for merge
#mn2: median of the second segment to test for merge
#pv.thres: combine two segments if p-value for wilcoxon rank sum test between
#          the observed values falling into two segments is greater than  pv.thres
#OR
#thresAbs: combine two segments if difference between their medians (copy number
#           transformed) is less than or equal to thresAbs

combine.func=function(diff,vecObs, vecPredNow, mnNow, mn1, mn2, pv.thres=0.0001, thresAbs=0)
{
  #observed values in the first segment
  vec1=vecObs[which(vecPredNow==mn1)]
  #observed values in the second segment
  vec2=vecObs[which(vecPredNow==mn2)]
  
  #if difference between segment medians does not exceed thresAbs, then set pv=1
  if (diff<=thresAbs) {
    pv=1
  }
  #otherwise test for difference in mean based on observed values
  else {
    if((length(vec1) > 10 & length(vec2) > 10) | sum(length(vec1),length(vec2))>100){
      pv=wilcox.test(vec1,vec2)$p.value
    }
    else{pv=wilcox.test(vec1,vec2,exact=T)$p.value  }	#/10^max(mn1,mn2)
    if(length(vec1) <= 3 | length(vec2) <= 3){pv=0}		
  }
  index.merged<-numeric()
  #if p-value exceeds pv.thres
  if (pv > pv.thres) 	{
    #combine observed values
    vec=c(vec1,vec2)
    # Index values to be updated
    index.merged=which((vecPredNow==mn1) | (vecPredNow==mn2))		
    #update predicted values by median of the observed values
    vecPredNow[index.merged]=median(vec, na.rm=TRUE)
    #update segment medians  median of the observed values and remove one of the duplicates
    mnNow[which((mnNow==mn1) | (mnNow==mn2))]=median(vec, na.rm=TRUE)
    mnNow=unique(mnNow)
  }
  list(mnNow=mnNow, vecPredNow=vecPredNow, pv=pv)
}



########################

# pv.thres: threshold for wilcoxon test
# ansari.sign: significance threshold for ansari-bradley test

# Level Merging Function
MergeLevels <- function(vecObs,vecPred,pv.thres=0.0001,ansari.sign=0.05){
  
  # Initializing threshold vector for keeping track of thresholds
  sq<-numeric()
  
  #initializing threshold index (threshold count)
  j=0
  
  #initializing ansari p-values to keep track of ansari p-values for each threshold in sq
  ansari=numeric()
  
  # Initialize levels count
  lv=numeric()
  
  # Start with threshold 0.05, and flag=0 indicating significance not yet reached, backtracking not begun
  thresAbs=0.05
  flag=0
  while (1){
    j=j+1
    # Save current threshold
    sq[j]<-thresAbs
    
    # temporary predicted values (to be updated)
    vecPredNow=vecPred
    
    #unmissing unique segment medians
    mnNow=unique(vecPred)
    mnNow=mnNow[!is.na(mnNow)]
    
    #continuing indicator otherwise get out of the loop
    cont=0
    
    while(cont==0 & length(mnNow)>1) {
      
      mnNow=sort(mnNow)  #currennt sorted vector of means
      n <- length(mnNow)  # number of means in mnNow
      # cat("\r",n,":",length(unique(vecPred)))
      # Get distances translated to copy number differences
      # Only distances to closest levels
      d<-(2*2^mnNow)[-n]-(2*2^mnNow)[-1]
      
      #order distance between means with the closest on top and corresponding indices
      dst<-cbind(abs(d)[order(abs(d))],(2:n)[order(abs(d))],(1:(n-1))[order(abs(d))])
      
      #for each pair of means
      for (i in 1:nrow(dst)) 	{
        #set continuity index to "NOT continue" (=1)
        cont=1
        #test for combining of the two segment means
        out=combine.func(diff=dst[i,1],vecObs, vecPredNow, mnNow, mn1=mnNow[dst[i,2]], mn2=mnNow[dst[i,3]], pv.thres=pv.thres, thresAbs=thresAbs)
        #if combine?
        if (out$pv > pv.thres) {
          
          #set continuity index to "YES" (=0) and break out of the current pairs loop
          cont=0
          
          #update predicted values and segments
          vecPredNow=out$vecPredNow
          mnNow=out$mnNow
          break
        }		
      }		
    }
    # When done merging for a given threshold, test for significance
    ansari[j]=ansari.test(sort(vecObs-vecPredNow), sort(vecObs-vecPred))$p.value
    if(is.na(ansari[j])){ansari[j]=0}
    lv[j]=length(mnNow) # get number of levels
    
    # If p.value is less than the significance threshold, set backtracking flag=1 (backtracking on)
    if(ansari[j]<ansari.sign){
      flag=1
    }
    if(flag==2){ break } 
    
    # If backtracking is on, a smaller threshold is attempted
    if (flag){
      
      # If backtracking is on and p.value is higher than sign threshold, stop
      if (ansari[j]>ansari.sign | thresAbs == 0){
        
        # Don't merge at all if all tested threshold including 0 is significant
        if (ansari[j] <= ansari.sign) {
          vecPredNow=vecPred
          mnNow=unique(vecPred)
          mnNow=mnNow[!is.na(mnNow)]
        }
        break
      }
      
      # Attempt smaller threshold
      else {thresAbs=signif(thresAbs-0.005,3) }
    }
    else {thresAbs=thresAbs+0.1} # Increase threshold if backtracking is not on
    
    # Control step so function won't keep running, max threshold = 1 and if sign not reached, threshold = 0
    if (thresAbs >= 1){
      thresAbs=0
      flag=2
    }
  }
  
  
  # Return list of results
  return(list(vecMerged=vecPredNow,mnNow=mnNow,sq=sq,ansari=ansari))
}

# Apply merge levels to population segmentation data

apply_merge_levels <- function(bin_data,
                               popseg_long,
                               mccores = cpus) {
  
  message("Running merge levels.")
  
  long_data_t <- as.data.frame(t(popseg_long))
  
  # removing chromosomes info to calculate the ratios
  bin_data_cells <-
    bin_data %>% dplyr::select(-c(abspos, chrom, chrompos))
  
  bin_data_info <- 
    bin_data %>% dplyr::select(c(chrom, chrompos, abspos))
    
  ratio_data <-
    sweep(bin_data_cells, 2, apply(bin_data_cells, 2, mean), '/')
  names(ratio_data) <- names(long_data_t)

    merge_levels_data <- mclapply(seq_along(ratio_data), function(x) {
      this_ratio <- ratio_data[, x]
      this_seg <- long_data_t[, x]
      merge.obj <-
        MergeLevels(log(this_ratio + 1e-3, base = 2), this_seg)
      return(merge.obj$vecMerged)
    }, mc.cores = mccores)
    
  names(merge_levels_data) <- names(ratio_data)
  
  merge_levels_df <- bind_rows(merge_levels_data, .id = "cell")
  
  ratio_data_l <- ratio_data
  
  long_data_t <- as.data.frame(2 ^ merge_levels_df)
  
  long_data_t <- bind_cols(bin_data_info, long_data_t)
  
  write.table(long_data_t,
              file = o,
              row.names = F,
              sep = "\t")
  
}


# Capturing args and checks

args <- parser$parse_args()

bin <- args$bincounts
o <- args$output
popseg <- args$popseg
cpus <- args$cores


if (is.null(o)) stop("No output filename provided, use flag -o")

if (file.access(bin) == -1) {
  stop(sprintf("File (%s) does not exist.", bin))
} else {
  bin <- bin
}

if (file.access(popseg) == -1) {
  stop(sprintf("File (%s) does not exist.", popseg))
} else {
  popseg <- popseg
}


# Running merge levels

message("Reading data.")

inbin <- read_delim(bin,
                    delim = "\t")

inpopseg <- read_delim(popseg,
                       delim = "\t")

popseg_long <- create_popseg_long(inpopseg)

apply_merge_levels(bin_data = inbin,
                   popseg_long = popseg_long,
                   mccores = cpus)

message("Done.")
