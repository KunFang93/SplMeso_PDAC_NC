# The splanchnic mesenchyme is the tissue of origin for fibroblasts in the pancreas during homeostasis and tumorigenesis
Kun Part of scRNA-seq analysis codes for the paper

## Data
Seurat oject used in codes should be found under GSE200903, [FibroXplorer](https://www.fibroxplorer.com) and [Tabula Sapiens](https://figshare.com/projects/Tabula_Sapiens/100973)

## System requirements
R with packages
```
library(Seurat)
library(SeuratDisk)
library(sctransform)
library(patchwork)
library(dplyr)
library(ggplot2)
library(reshape)
library(writexl)
library(readxl)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(collections)
library(Matrix)
```
**SessionInfo**
```
R version 4.2.1 (2022-06-23)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Monterey 12.5

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] Matrix_1.4-1          collections_0.3.5     circlize_0.4.15       forcats_0.5.1        
 [5] stringr_1.4.0         purrr_0.3.4           readr_2.1.2           tidyr_1.2.0          
 [9] tibble_3.1.8          tidyverse_1.3.2       ComplexHeatmap_2.12.0 readxl_1.4.0         
[13] writexl_1.4.0         reshape_0.8.9         ggplot2_3.3.6         dplyr_1.0.9          
[17] patchwork_1.1.1       sctransform_0.3.3     SeuratDisk_0.0.0.9020 sp_1.5-0             
[21] SeuratObject_4.1.0    Seurat_4.1.1         

loaded via a namespace (and not attached):
  [1] backports_1.4.1       plyr_1.8.7            igraph_1.3.4          lazyeval_0.2.2       
  [5] splines_4.2.1         listenv_0.8.0         scattermore_0.8       digest_0.6.29        
  [9] foreach_1.5.2         htmltools_0.5.3       fansi_1.0.3           magrittr_2.0.3       
 [13] tensor_1.5            googlesheets4_1.0.0   cluster_2.1.3         doParallel_1.0.17    
 [17] ROCR_1.0-11           tzdb_0.3.0            globals_0.15.1        modelr_0.1.8         
 [21] matrixStats_0.62.0    spatstat.sparse_2.1-1 colorspace_2.0-3      rvest_1.0.2          
 [25] ggrepel_0.9.1         haven_2.5.0           crayon_1.5.1          jsonlite_1.8.0       
 [29] progressr_0.10.1      spatstat.data_2.2-0   survival_3.3-1        zoo_1.8-10           
 [33] iterators_1.0.14      glue_1.6.2            polyclip_1.10-0       gtable_0.3.0         
 [37] gargle_1.2.0          leiden_0.4.2          GetoptLong_1.0.5      future.apply_1.9.0   
 [41] shape_1.4.6           BiocGenerics_0.42.0   abind_1.4-5           scales_1.2.0         
 [45] DBI_1.1.3             spatstat.random_2.2-0 miniUI_0.1.1.1        Rcpp_1.0.9           
 [49] viridisLite_0.4.0     xtable_1.8-4          clue_0.3-61           reticulate_1.25      
 [53] spatstat.core_2.4-4   bit_4.0.4             stats4_4.2.1          htmlwidgets_1.5.4    
 [57] httr_1.4.3            RColorBrewer_1.1-3    ellipsis_0.3.2        ica_1.0-3            
 [61] pkgconfig_2.0.3       uwot_0.1.11           dbplyr_2.2.1          deldir_1.0-6         
 [65] utf8_1.2.2            tidyselect_1.1.2      rlang_1.0.4           reshape2_1.4.4       
 [69] later_1.3.0           munsell_0.5.0         cellranger_1.1.0      tools_4.2.1          
 [73] cli_3.3.0             generics_0.1.3        broom_1.0.0           ggridges_0.5.3       
 [77] fastmap_1.1.0         goftest_1.2-3         bit64_4.0.5           fs_1.5.2             
 [81] fitdistrplus_1.1-8    RANN_2.6.1            pbapply_1.5-0         future_1.27.0        
 [85] nlme_3.1-158          mime_0.12             xml2_1.3.3            hdf5r_1.3.5          
 [89] compiler_4.2.1        rstudioapi_0.13       plotly_4.10.0         png_0.1-7            
 [93] spatstat.utils_2.3-1  reprex_2.0.1          stringi_1.7.8         rgeos_0.5-9          
 [97] lattice_0.20-45       vctrs_0.4.1           pillar_1.8.0          lifecycle_1.0.1      
[101] spatstat.geom_2.4-0   lmtest_0.9-40         GlobalOptions_0.1.2   RcppAnnoy_0.0.19     
[105] data.table_1.14.2     cowplot_1.1.1         irlba_2.3.5           httpuv_1.6.5         
[109] R6_2.5.1              promises_1.2.0.1      KernSmooth_2.23-20    gridExtra_2.3        
[113] IRanges_2.30.0        parallelly_1.32.1     codetools_0.2-18      MASS_7.3-58          
[117] assertthat_0.2.1      rjson_0.2.21          withr_2.5.0           S4Vectors_0.34.0     
[121] hms_1.1.1             mgcv_1.8-40           parallel_4.2.1        rpart_4.1.16         
[125] googledrive_2.0.0     Rtsne_0.16            lubridate_1.8.0       shiny_1.7.2          

```
## Installation guide
All dependent packages can be install with [bioconductor](https://www.bioconductor.org) or install.packages() command.  
Installation time varies based on the computer. Should not be longer than 30 mins. 


## Demo
**human_correlation.R** will generate heatmap of genes highly expressed in both the splanchnic mesenchyme and human pancreatic CAFs compared to human 
pancreatic TRFs (upper block), and genes highly expressed in both the splanchnic mesenchyme and human pancreatic TRFs 
compared to human pancreatic CAFs (lower block). Spl, splanchnic; Meso, mesoderm; CAFs, cancer associated fibroblasts; 
TRFs, tissue resident fibroblasts.

<img src="https://github.com/KunFang93/SplMeso_PDAC_NC/blob/main/Spl_Meos_human.PDAC.Normal_072922.png" width="900">

**mouse_analysis.R** will generate 1. the heatmap of genes highly expressed in both the splanchnic mesenchyme and CAFs compared 
to TRFs (upper block), and genes highly expressed in both the splanchnic mesenchyme and TRFs compared to CAFs (lower block). 2. generate Dimplot for Eng.
3. Find DEGs of 6 vs 0,1,13 clusters  and 4.UMAPs for different subtypes of fibroblasts in between IcreT and KPFIcreT tissues 

*1*.
<img src="https://github.com/KunFang93/SplMeso_PDAC_NC/blob/main/twogroup_FCordered_072222.png" width="900">


*2*.
<img src="https://github.com/KunFang93/SplMeso_PDAC_NC/blob/main/Eng_dimplot.png" width="900">

*4*.
<img src="https://github.com/KunFang93/SplMeso_PDAC_NC/blob/main/UMAPs.png" width="900">


## Instructions for use
It is recommeded to run codes in Rstudio.
Runing the code line by line :)
All pics should be reproducible from the codes


## Contact Us
Kun Fang: fangk@mcw.edu; Lu Han: hanl@musc.edu

