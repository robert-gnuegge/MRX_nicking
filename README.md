# MRX_nicking [![DOI](https://zenodo.org/badge/438784647.svg)](https://zenodo.org/badge/latestdoi/438784647)

This repository contains `.bash` and `.R` scripts used in our MRX nicking paper.

See publication for command line software versions. Find R package versions below (output of `sessionInfo()`):
```
R version 4.2.2 Patched (2022-11-10 r83330)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 18.04.6 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1

locale:
 [1] LC_CTYPE=en_US.UTF-8          LC_NUMERIC=C                  LC_TIME=en_US.UTF-8           LC_COLLATE=en_US.UTF-8        LC_MONETARY=en_US.UTF-8      
 [6] LC_MESSAGES=en_US.UTF-8       LC_PAPER=en_US.UTF-8          LC_NAME=en_US.UTF-8           LC_ADDRESS=en_US.UTF-8        LC_TELEPHONE=en_US.UTF-8     
[11] LC_MEASUREMENT=en_US.UTF-8    LC_IDENTIFICATION=en_US.UTF-8

attached base packages:
 [1] parallel  grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] DNAshapeR_1.24.0                        flowViz_1.60.2                          lattice_0.20-45                         flowCore_2.8.0                         
 [5] BiocManager_1.30.19                     beanplot_1.3.1                          propagate_1.0-6                         minpack.lm_1.2-2                       
 [9] ff_4.0.7                                bit_4.0.5                               Rcpp_1.0.9                              tmvtnorm_1.5                           
[13] gmm_1.7                                 sandwich_3.0-2                          Matrix_1.5-3                            mvtnorm_1.1-3                          
[17] MASS_7.3-58.1                           tictoc_1.1                              rmelting_1.12.0                         plotrix_3.8-2                          
[21] BSgenome.Scerevisiae.UCSC.sacCer3_1.4.0 BSgenome_1.64.0                         Gviz_1.40.1                             GenomicAlignments_1.32.0               
[25] Rsamtools_2.12.0                        SummarizedExperiment_1.26.1             Biobase_2.56.0                          MatrixGenerics_1.8.0                   
[29] matrixStats_0.63.0                      rtracklayer_1.56.1                      Biostrings_2.64.0                       XVector_0.36.0                         
[33] GenomicRanges_1.48.0                    GenomeInfoDb_1.32.1                     IRanges_2.30.0                          S4Vectors_0.34.0                       
[37] BiocGenerics_0.42.0                    

loaded via a namespace (and not attached):
  [1] spam_2.9-1               backports_1.4.1          Hmisc_4.7-2              BiocFileCache_2.4.0      lazyeval_0.2.2           splines_4.2.2           
  [7] BiocParallel_1.30.0      ggplot2_3.4.0            digest_0.6.30            ensembldb_2.20.1         htmltools_0.5.3          viridis_0.6.2           
 [13] fansi_1.0.3              magrittr_2.0.3           checkmate_2.1.0          memoise_2.0.1            cluster_2.1.4            RcppParallel_5.1.5      
 [19] cytolib_2.8.0            prettyunits_1.1.1        jpeg_0.1-10              colorspace_2.0-3         blob_1.2.3               rappdirs_0.3.3          
 [25] rbibutils_2.2.10         xfun_0.35                dplyr_1.0.10             crayon_1.5.2             RCurl_1.98-1.9           hexbin_1.28.2           
 [31] survival_3.4-0           VariantAnnotation_1.42.0 zoo_1.8-11               glue_1.6.2               gtable_0.3.1             zlibbioc_1.42.0         
 [37] DelayedArray_0.22.0      IDPmisc_1.1.20           maps_3.4.1               scales_1.2.1             DBI_1.1.3                viridisLite_0.4.1       
 [43] progress_1.2.2           htmlTable_2.4.1          foreign_0.8-83           dotCall64_1.0-2          Formula_1.2-4            htmlwidgets_1.5.4       
 [49] httr_1.4.4               RColorBrewer_1.1-3       ellipsis_0.3.2           pkgconfig_2.0.3          XML_3.99-0.13            rJava_1.0-6             
 [55] nnet_7.3-18              dbplyr_2.2.1             deldir_1.0-6             utf8_1.2.2               tidyselect_1.2.0         rlang_1.0.6             
 [61] AnnotationDbi_1.58.0     munsell_0.5.0            tools_4.2.2              cachem_1.0.6             cli_3.4.1                generics_0.1.3          
 [67] RSQLite_2.2.19           stringr_1.5.0            fastmap_1.1.0            yaml_2.3.6               knitr_1.41               bit64_4.0.5             
 [73] KEGGREST_1.36.0          AnnotationFilter_1.20.0  xml2_1.3.3               biomaRt_2.52.0           compiler_4.2.2           rstudioapi_0.14         
 [79] filelock_1.0.2           curl_4.3.3               png_0.1-8                tibble_3.1.8             stringi_1.7.8            GenomicFeatures_1.48.0  
 [85] fields_14.1              ProtGenerics_1.28.0      vctrs_0.5.1              pillar_1.8.1             lifecycle_1.0.3          Rdpack_2.4              
 [91] data.table_1.14.6        bitops_1.0-7             R6_2.5.1                 BiocIO_1.6.0             latticeExtra_0.6-30      KernSmooth_2.23-20      
 [97] gridExtra_2.3            RProtoBufLib_2.8.0       dichromat_2.0-0.1        assertthat_0.2.1         rjson_0.2.21             GenomeInfoDbData_1.2.8  
[103] hms_1.1.2                rpart_4.1.19             biovizBase_1.44.0        base64enc_0.1-3          interp_1.1-3             restfulr_0.0.15  
``` 

This project is licensed under the terms of the MIT license.
