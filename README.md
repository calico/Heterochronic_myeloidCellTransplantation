# Heterochronic Brain Myeloid Cell Transplantation analysis
This GitHub repository documents the core analyses of scRNA-seq and Spatial-seq data of microglia and heterochronic myeloid cell transplantation experiments. 

## Background
Aging, the key risk factor for cognitive decline, impacts the brain in a region- and cell type-specific manner. Microglia are considered among the fastest aging cell types; however, it remains unclear whether this is intrinsically mediated or is driven by age-related changes in neighboring cells. Here, we describe analyses of scRNA-seq and Spatial-seq data derived from in vivo heterochronic myeloid cell replacement experiments. 

<img width="3182" height="1081" alt="Fig_Shiny" src="https://github.com/user-attachments/assets/75055703-fa22-430a-8a3d-4954d1ac5ecf" />


We document the following:

* Loading, quality control, integration and clustering of scRNA-seq or Spatial-seq (CosMx) data
* Annotation and filtering of cell clusters (scRNA-seq)
* Regional clustering and annotation of cells (Spatial-seq)
* Differential gene expression analysis across regions and/or age groups
* Definition of a set of differentially expressed genes that shift during aging and under heterochronic reconstitution
* Perform differential cell-cell interaction analysis
* Co-integration Spatial- and scRNA-seq data
* Cell-cell neighborhood analysis in Spatial-seq data
* Define gene signatures and quantify signature scores in Spatial- and scRNA-seq data


## Associated Manuscript
If using data or scripts of this study, please cite the following pre-print:
Heterochronic myeloid cell replacement reveals the local brain environment as key driver of microglia aging.
Claire Gizowski, Galina Popova, Heather Shin, Marius M Mader, Wendy Craft, Bernd J Wranik, Wenjun Kong, Yuheng C Fu, Constanze Depp, Tzuhua D Lin, Baby Martin-McNulty, Han Tai, Nicole Fong, Devyani Jogran, Kayla Leung, Agnieszka Wendorff, David Hendrickson, Astrid Gillich, Andy Chang, Beth Stevens, Marius Wernig, Oliver Hahn. bioRxiv. doi: XXX

## Data availability
The sequencing datasets analyzed during the current study are available in the BioProject repository under accession number PRJNA600501, and in the Gene Expression Omnibus repository under accession numbers GSE306331, GSE306335, GSE306344, GSE306351 and GSE306354.


# Installation Instructions
Please install the R packages listed in the following 'Dependencies' section.

## Dependencies
The original analyses were peformed within the following enviornment: 

R 4.4.0 (2024-04-24) Platform: x86_64-pc-linux-gnu Running under: Ubuntu 24.04.2 LTS

The following R packages are required for this script:

* anndata
* CellChat
* data.table
* dbscan
* dplyr
* future
* GEOquery
* ggplot2
* ggrepel
* harmony
* limma
* lsmeans
* Matrix
* parallel
* patchwork
* plyr
* purrr
* RANN
* readr
* scales
* Seurat
* stringr
* VISION
