# Ontogeny-independent expression of LPCAT2 in granuloma macrophages during experimental visceral leishmaniasis

## Authors
Shoumit Dey<sup>1#</sup>, Jian-Hua Cao<sup>2#</sup>, Benjamin Balluff<sup>2</sup>, Nidhi Sharma Dey<sup>1</sup>, Lesley Gilbert<sup>3</sup>, Sally James<sup>3</sup>, Adam A. Dowle<sup>3</sup>, Grant Calder<sup>3</sup>, Peter O'Toole<sup>3</sup>, Ron M. A. Heeren<sup>2*</sup>, Paul M. Kaye<sup>1*</sup>

<sup>1</sup> York Biomedical Research Institute, Hull York Medical School, University of York, UK  
<sup>2</sup> Maastricht MultiModal Molecular Imaging (M4I) Institute, Maastricht University, the Netherlands  
<sup>3</sup> Biosciences Technology Facility, Department of Biology, University of York, UK  
<sup>#</sup> These authors contributed equally.  
<sup>*</sup> Correspondence: paul.kaye@york.ac.uk; r.heeren@maastrichtuniversity.nl

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

### Key Python packages required:
- scanpy
- anndata
- cell2location
- scvi-tools
- pandas
- numpy
- matplotlib
- seaborn

A `requirements.txt` for Python dependencies is provided in this repository.

<details>
<summary>sessionInfo() (click to expand)</summary>

```
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
```

</details>

### Running Time:
- End-to-end analysis (R scripts): ~1–2 hours on Windows 10 with 64GB RAM, 8-core CPU (Intel Core i7-9700 or equivalent).

### Additional files required:
Rds files will be made available upon request.

## Repository Structure & Instructions

### File overview

| File | Description |
|------|-------------|
| `scRNAseq_integrated_prepare.Rmd` | scRNA-seq data loading, QC, integration, and clustering |
| `scRNAseq_downstream_analysis.Rmd` | scRNA-seq downstream analysis, DE, and figure generation |
| `spatial_integrated_prepare_RNA.Rmd` | Spatial transcriptomics (Visium) data loading, QC, integration, and clustering |
| `spatial_integrated_prepare_MSI.Rmd` | Mass spectrometry imaging (MSI) data loading, integration, and clustering |
| `spatial_downstream_analysis.Rmd` | Spatial multi-modal downstream analysis and figure generation |
| `config.py` | Cell2location pipeline configuration and entry point |
| `train_ref.py` | Cell2location reference model training |
| `load_ref.py` | Load trained reference model and export signatures |
| `load_query.py` | Load and QC spatial (Visium) query data |
| `deconvolute.py` | Cell2location spatial mapping and deconvolution |
| `nmf_compartments.py` | NMF-based co-location analysis of cell types |
| `sd22.5.1.sh` | SLURM submission script for cell2location pipeline |
| `sd22.5.1.nmf.sh` | SLURM submission script for NMF co-location analysis |
| `utils.R` | Shared R utility functions (valley finding, correlation) |

### Single-cell RNA-seq analysis:
- To reproduce clustering and annotation, start with `scRNAseq_integrated_prepare.Rmd`.
- Alternatively, use prepared .rds files (available on Zenodo) and start directly with `scRNAseq_downstream_analysis.Rmd`.

### Spatial Transcriptomics (Visium):
- To reproduce spatial integration, start with `spatial_integrated_prepare_RNA.Rmd` (for transcriptomics) and `spatial_integrated_prepare_MSI.Rmd` (for mass spectrometry imaging).
- Alternatively, load prepared .rds files and start directly with `spatial_downstream_analysis.Rmd`.

### Cell2location Spatial Integration:
- Cell2location was run on GPU nodes using the University of York's HPC (Viking), requiring ~2 hours for complete spatial mapping.
- Python scripts provided (`config.py`, `train_ref.py`, `load_ref.py`, `load_query.py`, `deconvolute.py`) with necessary .h5ad files available upon publication on Zenodo.
- NMF co-location analysis is run separately via `nmf_compartments.py`.
- SLURM job submission scripts (`sd22.5.1.sh`, `sd22.5.1.nmf.sh`) are provided as templates for HPC usage.

## Quickstart (Example workflow)

### Option 1: Raw data analysis (post-publication)
- Download raw Visium and scRNA-seq data from GEO accession GSE290324 and GSE290325 respectively.
- Run analysis notebooks sequentially from data preparation through figure generation.

### Option 2: Prepared data analysis
- Download .rds files and necessary CSV metadata from Zenodo.
- Begin directly from `scRNAseq_downstream_analysis.Rmd` and `spatial_downstream_analysis.Rmd`.

## Contact and Correspondence
- Shoumit Dey : shoumit.dey@york.ac.uk / shoumit@gmail.com

## License
This project is covered under the MIT License.
