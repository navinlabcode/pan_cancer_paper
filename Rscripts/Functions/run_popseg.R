#' Population Segmentation from VarBin bin count inputs with multipcf
#'
#' @param inuber A bin count matrix resulting from the VarBin pipeline
#' @param outprefix Output file name
#' @param ingamma Gamma penalty, see multipcf help, default is 20.
#' @param inpseudo A pseudonumber that can be added to the bin counts
#'
#' @author Alexander Davis, Darlan Conterno Minussi, Hanghui Ye
#' @return
#' @export
#' @examples
#' 
#'

# Setup
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(copynumber))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(reshape2))

parser <- ArgumentParser()

parser$add_argument("-i", "--inuber",
                    help = "A bin count matrix resulting from the VarBin pipeline")

parser$add_argument("-o", "--outprefix", type = "character",
                    help = "Output file name")

parser$add_argument("-g", "--ingamma", type = "integer",
                    help = "Gamma penalty, see multipcf help, default is 20.",
                    default = 20)

parser$add_argument("-p", "--inpseudo", type = "integer",
                    help = "A pseudonumber that can be added to the bin counts.",
                    default = 0)

# Population segmentation function

run_popseg <- function(inuber, outprefix, ingamma = 20, inpseudo = 0) {
  
  # Load in the bincount uber and turn it into logratios with chrom and chrompos columns
  uber.matrix <- 
    read_tsv(inuber,
             col_types=cols(chrom=col_character(), chrompos=col_double(), abspos=col_double(),
                            .default=col_double()))
  
  # checks
  if (any(uber.matrix < 0) == TRUE) {
    stop("Population segmenter needs a bin count matrix, your input contains negative values, please make sure you have the correct matrix")
  }
  
  if (length(grep("chrompos", names(uber.matrix))) > 1) {
    stop("You have duplicated chrompos columns. Execution stopped")
  }
  
  if (length(grep("abspos", names(uber.matrix))) > 1) {
    stop("You have duplicated abspos columns. Execution stopped")
  }
  
  bincount.mat <- as.matrix(dplyr::select(uber.matrix, -chrom, -chrompos, -abspos))+inpseudo
  ratio.mat <- sweep(bincount.mat, 2, apply(bincount.mat, 2, mean), '/')
  annotated.log.gc.ratio <-
    bind_cols(dplyr::select(uber.matrix, chrom, chrompos),
              as.data.frame(log2(ratio.mat)))
  
  # Winsorize (if the data is beyond some cutoff compared to its neighbors, set it
  # equal to that cutoff)
  winsorized.logratios <- winsorize(as.data.frame(annotated.log.gc.ratio))
  # Segment the data
  # This is another place you might want to customize: different values of gamma
  # give you different segmentations
  seg.output <-
    multipcf(winsorized.logratios, gamma=ingamma,
             return.est=FALSE, assembly='hg19', fast=TRUE)
  
  readr::write_tsv(seg.output, sprintf("%s", outprefix))
  
}


# running population segmentation

message("Reading data.")

args <- parser$parse_args()

m <- args$inuber
o <- args$outprefix
g <- args$ingamma
p <- args$inpseudo

if (is.null(o)) stop("No output filename provided, use flag -o")

if (file.access(m) == -1) {
  stop(sprintf("File (%s) does not exist.", m))
} else {
  inuber <- m
}

run_popseg(inuber = inuber,
           outprefix = o,
           ingamma = g,
           inpseudo = p)

message("Done.")

