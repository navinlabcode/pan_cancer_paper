# A pan-cancer single-cell analysis of intratumoral copy number diversity and evolution 

## Instructions

This repository contains the R scripts used to create the main figures for the manuscript titled: **A pan-cancer single-cell analysis of intratumoral copy number diversity and evolution**. 

The `Rscripts` folder is divided into two subfolders: `Functions` and `Figures`. `Functions` includes the general functions used for the CNA analysis in the paper. `Figures` contains detailed step-by-step R code to generate the plots for main figures 1-6. The `Data` folder includes the processed datasets used to plot these figures.

The raw sequencing data is available at the Sequence Read Archive (SRA) under the accession number PRJNA1013415.

## Dependencies

Part of the R scripts were adapted from [*CopyKit*](https://github.com/navinlabcode/copykit) (v0.1.0). *CopyKit* can be installed using the following command:

```r
devtools::install_github(repo = "navinlabcode/copykit", ref = "f709a48")
```
Session info:
```
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux

Matrix products: default
BLAS/LAPACK: /usr/lib64/libopenblasp-r0.3.3.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] grid      stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] VennDiagram_1.7.1                       futile.logger_1.4.3                     ggtree_3.2.1                           
 [4] ComplexHeatmap_2.10.0                   Homo.sapiens_1.3.1                      TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
 [7] org.Hs.eg.db_3.14.0                     GO.db_3.14.0                            OrganismDbi_1.36.0                     
[10] GenomicFeatures_1.46.5                  AnnotationDbi_1.56.2                    Biobase_2.54.0                         
[13] GenomicRanges_1.46.1                    GenomeInfoDb_1.30.1                     IRanges_2.28.0                         
[16] S4Vectors_0.32.3                        BiocGenerics_0.40.0                     maftools_2.10.05                       
[19] ape_5.6-2                               scatterpie_0.1.8                        stringr_1.5.0                          
[22] survival_3.3-0                          survminer_0.4.9                         ggpubr_0.4.0                           
[25] readxl_1.3.1                            dplyr_1.0.8                             ggplot2_3.3.6.9000                     
[28] reticulate_1.28                        

loaded via a namespace (and not attached):
  [1] utf8_1.2.3                  R.utils_2.11.0              tidyselect_1.2.0            RSQLite_2.2.10             
  [5] BiocParallel_1.28.3         munsell_0.5.0               codetools_0.2-18            withr_2.5.0                
  [9] colorspace_2.1-0            filelock_1.0.2              knitr_1.40                  rstudioapi_0.13            
 [13] ggsignif_0.6.3              labeling_0.4.2              MatrixGenerics_1.6.0        GenomeInfoDbData_1.2.7     
 [17] KMsurv_0.1-5                polyclip_1.10-4             bit64_4.0.5                 farver_2.1.1               
 [21] vctrs_0.5.1                 treeio_1.18.1               generics_0.1.3              lambda.r_1.2.4             
 [25] xfun_0.37                   BiocFileCache_2.2.1         R6_2.5.1                    doParallel_1.0.17          
 [29] clue_0.3-60                 pals_1.7                    bitops_1.0-7                cachem_1.0.6               
 [33] gridGraphics_0.5-1          DelayedArray_0.20.0         assertthat_0.2.1            BiocIO_1.4.0               
 [37] scales_1.2.0.9000           gtable_0.3.1                rlang_1.0.6                 GlobalOptions_0.1.2        
 [41] splines_4.1.2               rtracklayer_1.54.0          rstatix_0.7.0               lazyeval_0.2.2             
 [45] dichromat_2.0-0             prismatic_1.1.0             broom_0.7.12                BiocManager_1.30.18        
 [49] yaml_2.3.7                  abind_1.4-5                 backports_1.4.1             RBGL_1.70.0                
 [53] tools_4.1.2                 ggplotify_0.1.0             ellipsis_0.3.2              RColorBrewer_1.1-3         
 [57] Rcpp_1.0.10.2               progress_1.2.2              zlibbioc_1.40.0             purrr_0.3.4                
 [61] RCurl_1.98-1.6              prettyunits_1.1.1           GetoptLong_1.0.5            cowplot_1.1.1              
 [65] zoo_1.8-11                  SummarizedExperiment_1.24.0 haven_2.4.3                 cluster_2.1.2              
 [69] magrittr_2.0.3              magick_2.7.3                data.table_1.14.2           futile.options_1.0.1       
 [73] openxlsx_4.2.5              circlize_0.4.14             amap_0.8-18                 matrixStats_0.63.0         
 [77] evaluate_0.20               hms_1.1.1                   patchwork_1.1.2             xtable_1.8-4               
 [81] XML_3.99-0.10               rio_0.5.29                  gridExtra_2.3               shape_1.4.6                
 [85] compiler_4.1.2              biomaRt_2.50.3              tibble_3.1.8                maps_3.4.0                 
 [89] crayon_1.5.2                R.oo_1.24.0                 htmltools_0.5.2             mgcv_1.8-39                
 [93] ggfun_0.0.6                 tidyr_1.2.0                 aplot_0.1.2                 lubridate_1.8.0            
 [97] DBI_1.1.2                   tweenr_1.0.2                formatR_1.11                dbplyr_2.1.1               
[101] MASS_7.3-55                 rappdirs_0.3.3              Matrix_1.5-3                car_3.0-11                 
[105] cli_3.6.0                   R.methodsS3_1.8.1           forcats_0.5.1               pkgconfig_2.0.3            
[109] km.ci_0.5-2                 GenomicAlignments_1.30.0    foreign_0.8-82              xml2_1.3.3                 
[113] paletteer_1.4.0             foreach_1.5.2               XVector_0.34.0              snakecase_0.11.0           
[117] yulab.utils_0.0.4           digest_0.6.30               janitor_2.1.0               graph_1.72.0               
[121] Biostrings_2.62.0           rmarkdown_2.20              cellranger_1.1.0            survMisc_0.5.5             
[125] tidytree_0.3.8              restfulr_0.0.15             curl_4.3.2                  Rsamtools_2.10.0           
[129] rjson_0.2.21                lifecycle_1.0.3             nlme_3.1-155                jsonlite_1.8.4             
[133] carData_3.0-5               mapproj_1.2.8               fansi_1.0.3                 pillar_1.8.1               
[137] lattice_0.20-45             KEGGREST_1.34.0             fastmap_1.1.0               httr_1.4.4                 
[141] glue_1.6.2                  zip_2.2.0                   png_0.1-8                   iterators_1.0.14           
[145] bit_4.0.4                   ggforce_0.3.3               stringi_1.7.12              rematch2_2.1.2             
[149] blob_1.2.2                  memoise_2.0.1  

```
## Contact

For any additional information, please contact the corresponding author via [email](mailto:nnavin@mdanderson.org). 
