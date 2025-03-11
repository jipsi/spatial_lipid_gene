# Ontogeny-independent expression of LPCAT2 in granuloma macrophages during experimental visceral leishmaniasis

## Authors
Shoumit Dey<sup>1#</sup>, Jian-Hua Cao<sup>2#</sup>, Benjamin Balluff<sup>2</sup>, Nidhi Sharma Dey<sup>1</sup>, Lesley Gilbert<sup>3</sup>, Sally James<sup>3</sup>, Adam A. Dowle<sup>3</sup>, Grant Calder<sup>3</sup>, Peter O'Toole<sup>3</sup>, Ron M. A. Heeren<sup>2*</sup>, Paul M. Kaye<sup>1*</sup>

<sup>1</sup> York Biomedical Research Institute, Hull York Medical School, University of York, UK  
<sup>2</sup> Maastricht MultiModal Molecular Imaging (M4I) Institute, Maastricht University, the Netherlands  
<sup>3</sup> Biosciences Technology Facility, Department of Biology, University of York, UK  
<sup>#</sup> These authors contributed equally.  
<sup>*</sup> Correspondence: paul.kaye@york.ac.uk; r.heeren@maastrichtuniversity.nl

Published in: Nature Communications

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
- Download raw scRNA-seq and Visium data from GEO accession GSEXXXXX.
- Run analysis notebooks sequentially from data preparation through figure generation.

### Option 2: Prepared data analysis
- Download .rds files and necessary CSV metadata from Zenodo.
- Begin directly from `scRNAseq_downstream_analysis.Rmd` and `spatial_downstream_analysis.Rmd`.

## Notes on Experimental Design
- Animal model: Female C57BL/6J mice were used exclusively to minimize variability related to sex-dependent immune responses and disease progression during visceral leishmaniasis.

## Contact and Correspondence
- Paul Kaye: paul.kaye@york.ac.uk
- Ron Heeren: r.heeren@maastrichtuniversity.nl
- Shoumit Dey (data/code inquiries): shoumit.dey@york.ac.uk

## License
This project is covered under the MIT License.
