# README #

### Analysing omics data sets with weighted nodes networks (WNNets)

###### G. Tosadori, D. Di Silvestre, F. Spoto, P. Mauri, C. Laudanna and G. Scardoni

The code in the repository allows reproducing all the findings and all the figures that are presented in the published manuscript.

All the code was developed for R.

It is recommended, to obtain all the plots, to run the command from the **R command line** with the following syntax:

`source("scriptName.R", print.eval=T)`

or from RStudio using **Source with Echo** or Shift+Ctrl+Enter.

This is `sessionInfo()` loaded from the `03_results_analysis.R` script.

```
R version 4.1.0 (2021-05-18)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=it_IT.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=it_IT.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=it_IT.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=it_IT.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
[1] extrafont_0.17     igraph_1.2.6       RColorBrewer_1.1-2 ggrepel_0.9.1     
[5] UpSetR_1.4.0       reshape2_1.4.4     ggplot2_3.3.3     

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.6       pillar_1.5.1     compiler_4.1.0   plyr_1.8.6      
 [5] tools_4.1.0      lifecycle_1.0.0  tibble_3.1.0     gtable_0.3.0    
 [9] pkgconfig_2.0.3  rlang_0.4.10     DBI_1.1.1        Rttf2pt1_1.3.8  
[13] gridExtra_2.3    withr_2.4.1      dplyr_1.0.5      stringr_1.4.0   
[17] generics_0.1.0   vctrs_0.3.6      tidyselect_1.1.0 glue_1.4.2      
[21] R6_2.5.0         fansi_0.4.2      purrr_0.3.4      extrafontdb_1.0 
[25] magrittr_2.0.1   scales_1.1.1     ellipsis_0.3.1   assertthat_0.2.1
[29] colorspace_2.0-0 utf8_1.2.1       stringi_1.5.3    munsell_0.5.0   
[33] crayon_1.4.1    
```

STRING database v11 was used to obtain protein-protein interactions information.

Experimental only interactions were considered.

Medium confidence, i.e. 400+, threshold was used since its STRING's default.
