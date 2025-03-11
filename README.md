# Ontogeny-independent expression of LPCAT2 in granuloma macrophages during experimental visceral leishmaniasis

## Authors
Shoumit Dey<sup>1#</sup>, Jian-Hua Cao<sup>2#</sup>, Benjamin Balluff<sup>2</sup>, Nidhi Sharma Dey<sup>1</sup>, Lesley Gilbert<sup>3</sup>, Sally James<sup>3</sup>, Adam A. Dowle<sup>3</sup>, Grant Calder<sup>3</sup>, Peter O'Toole<sup>3</sup>, Ron M. A. Heeren<sup>2*</sup>, Paul M. Kaye<sup>1*</sup>

<sup>1</sup> York Biomedical Research Institute, Hull York Medical School, University of York, UK  
<sup>2</sup> Maastricht MultiModal Molecular Imaging (M4I) Institute, Maastricht University, the Netherlands  
<sup>3</sup> Biosciences Technology Facility, Department of Biology, University of York, UK  
<sup>#</sup> These authors contributed equally.  
<sup>*</sup> Correspondence: paul.kaye@york.ac.uk; r.heeren@maastrichtuniversity.nl

Unpublished

## Summary

Granulomas are structured immune lesions that form during chronic infections, such as visceral leishmaniasis. This study employed spatial transcriptomics, mass spectrometry imaging (MSI), single-cell RNA sequencing, and proteomics to characterize hepatic granulomas induced by Leishmania donovani infection in mice. A novel role for lysophosphatidylcholine acyltransferase 2 (LPCAT2)-mediated phospholipid remodeling was identified within granuloma-associated macrophages. LPCAT2 expression was independent of macrophage ontogeny and correlated with macrophage activation. Our findings demonstrate significant lipid-driven immunometabolic changes associated with granulomatous inflammation.

## System Requirements

### OS
- Windows: Windows 10 x64 (recommended)
- Mac
- Linux (e.g., CentOS, Ubuntu)

### Software
- R (≥4.2.2)
- RStudio (optional; recommended ≥2022.02.3+492)
- Python (≥3.8 for spatial transcriptomics integration via cell2location)

### Key R packages required:
- Seurat (≥4.3.0)
- ggplot2
- dplyr
- patchwork
- corrplot
- reshape2

sessionInfo()
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

Matrix products: default

locale:
[1] LC_COLLATE=English_United Kingdom.utf8  LC_CTYPE=English_United Kingdom.utf8   
[3] LC_MONETARY=English_United Kingdom.utf8 LC_NUMERIC=C                           
[5] LC_TIME=English_United Kingdom.utf8    

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] tidyr_1.2.1                 reshape2_1.4.4              VennDiagram_1.7.3           futile.logger_1.4.3        
 [5] stringr_1.5.0               corrplot_0.92               gplots_3.1.3                spatstat_3.0-2             
 [9] spatstat.linnet_3.0-3       spatstat.model_3.0-2        rpart_4.1.19                spatstat.explore_3.0-5     
[13] nlme_3.1-161                spatstat.random_3.0-1       spatstat.geom_3.0-3         spatstat.data_3.0-0        
[17] EnhancedVolcano_1.16.0      ggrepel_0.9.2               readxl_1.4.3                plotly_4.10.1              
[21] data.table_1.14.6           ggpubr_0.5.0                RColorBrewer_1.1-3          mclust_6.0.0               
[25] SingleCellExperiment_1.20.0 SummarizedExperiment_1.28.0 Biobase_2.58.0              GenomicRanges_1.50.2       
[29] GenomeInfoDb_1.34.6         IRanges_2.32.0              S4Vectors_0.36.1            BiocGenerics_0.44.0        
[33] MatrixGenerics_1.10.0       matrixStats_0.63.0          patchwork_1.1.2             ggplot2_3.4.2              
[37] sqldf_0.4-11                RSQLite_2.2.20              gsubfn_0.7                  proto_1.0.0                
[41] dplyr_1.0.10                SeuratObject_4.1.3          Seurat_4.3.0               

loaded via a namespace (and not attached):
  [1] utf8_1.2.2             reticulate_1.27        tidyselect_1.2.0       htmlwidgets_1.6.1     
  [5] Rtsne_0.16             munsell_0.5.0          codetools_0.2-18       ica_1.0-3             
  [9] chron_2.3-58           future_1.30.0          miniUI_0.1.1.1         withr_2.5.0           
 [13] colorspace_2.0-3       progressr_0.13.0       knitr_1.41             rstudioapi_0.14       
 [17] ROCR_1.0-11            ggsignif_0.6.4         tensor_1.5             listenv_0.9.0         
 [21] labeling_0.4.2         GenomeInfoDbData_1.2.9 polyclip_1.10-4        bit64_4.0.5           
 [25] farver_2.1.1           parallelly_1.33.0      vctrs_0.5.1            generics_0.1.3        
 [29] lambda.r_1.2.4         xfun_0.36              R6_2.5.1               bitops_1.0-7          
 [33] spatstat.utils_3.0-1   cachem_1.0.6           DelayedArray_0.24.0    promises_1.2.0.1      
 [37] scales_1.2.1           gtable_0.3.1           globals_0.16.2         goftest_1.2-3         
 [41] rlang_1.1.1            splines_4.2.2          rstatix_0.7.1          lazyeval_0.2.2        
 [45] broom_1.0.4            yaml_2.3.6             abind_1.4-5            backports_1.4.1       
 [49] httpuv_1.6.7           tools_4.2.2            tcltk_4.2.2            ellipsis_0.3.2        
 [53] ggridges_0.5.4         Rcpp_1.0.9             plyr_1.8.8             zlibbioc_1.44.0       
 [57] purrr_1.0.0            RCurl_1.98-1.9         deldir_1.0-6           pbapply_1.6-0         
 [61] cowplot_1.1.1          zoo_1.8-11             cluster_2.1.4          magrittr_2.0.3        
 [65] futile.options_1.0.1   scattermore_0.8        lmtest_0.9-40          RANN_2.6.1            
 [69] fitdistrplus_1.1-8     mime_0.12              evaluate_0.19          xtable_1.8-4          
 [73] gridExtra_2.3          compiler_4.2.2         tibble_3.1.8           KernSmooth_2.23-20    
 [77] crayon_1.5.2           htmltools_0.5.4        mgcv_1.8-41            later_1.3.0           
 [81] DBI_1.1.3              formatR_1.13           MASS_7.3-58.1          Matrix_1.5-3          
 [85] car_3.1-1              cli_3.6.0              parallel_4.2.2         igraph_1.3.5          
 [89] pkgconfig_2.0.3        sp_1.5-1               spatstat.sparse_3.0-0  XVector_0.38.0        
 [93] digest_0.6.31          sctransform_0.3.5      RcppAnnoy_0.0.20       rmarkdown_2.19        
 [97] cellranger_1.1.0       leiden_0.4.3           uwot_0.1.14            shiny_1.7.4           
[101] gtools_3.9.4           lifecycle_1.0.3        jsonlite_1.8.4         carData_3.0-5         
[105] viridisLite_0.4.1      fansi_1.0.3            pillar_1.9.0           lattice_0.20-45       
[109] fastmap_1.1.0          httr_1.4.6             survival_3.5-0         glue_1.6.2            
[113] png_0.1-8              bit_4.0.5              stringi_1.7.8          blob_1.2.3            
[117] caTools_1.18.2         memoise_2.0.1          irlba_2.3.5.1          future.apply_1.10.0  

### Running Time:
- End-to-end analysis (R scripts): ~1–2 hours on Windows 10 with 64GB RAM, 8-core CPU (Intel Core i7-9700 or equivalent).

### Additional files required:
Available on Zenodo: https://zenodo.org/record/XXXXX

## Repository Structure & Instructions

### Single-cell RNA-seq analysis:
- To reproduce clustering and annotation, start with `scRNAseq_integrated_prepare_rds.Rmd`.
- Alternatively, use prepared .rds files (available on Zenodo) and start directly with `scRNAseq_downstream_analysis.Rmd`.

### Spatial Transcriptomics (Visium):
- To reproduce spatial integration, start with `spatial_integrated_prepare_rds.Rmd`.
- Alternatively, load prepared .rds files and start directly with `spatial_downstream_analysis.Rmd`.

### Mass Spectrometry Imaging (MSI):
- Raw MSI data processing steps, lipid annotations, and spatial mapping scripts are provided in `MSI_analysis.Rmd`.

### Cell sorting and proteomics:
- Flow cytometry sorting strategy and downstream proteomic analyses scripts provided in `proteomics_analysis.Rmd`.

### Cell2location Spatial Integration:
- Cell2location was run on GPU nodes using the University of York's HPC (Viking), requiring ~2 hours for complete spatial mapping.
- Python scripts provided (`config.py`, `train_ref.py`, `deconvolute.py`) with necessary .h5ad files available upon publication on Zenodo.

## Quickstart (Example workflow)

### Option 1: Raw data analysis (post-publication)
- Download raw Visium and scRNA seq data from GEO accession GSE290324 and GSE290325 respectively.
- Run analysis notebooks sequentially from data preparation through figure generation.

### Option 2: Prepared data analysis
- Download .rds files and necessary CSV metadata from Zenodo.
- Begin directly from `scRNAseq_downstream_analysis.Rmd` and `spatial_downstream_analysis.Rmd`.

## Contact and Correspondence
- Paul Kaye: paul.kaye@york.ac.uk
- Ron Heeren: r.heeren@maastrichtuniversity.nl
- Shoumit Dey (data/code inquiries): shoumit.dey@york.ac.uk

## License
This project is covered under the MIT License.
