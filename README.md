# Heterochronic myeloid cell replacement
This GitHub repository documents the core analyses of scRNA-seq and Spatial-seq data of microglia and heterochronic myeloid cell transplantation experiments. 

## Background
Aging, the key risk factor for cognitive decline, impacts the brain in a region- and cell type-specific manner. Microglia are considered among the fastest aging cell types; however, it remains unclear whether this is intrinsically mediated or is driven by age-related changes in neighboring cells. Here, we describe analyses of scRNA-seq and Spatial-seq data derived from in vivo heterochronic myeloid cell replacement experiments. We document the following:

* Loading, quality control, integration and clustering of scRNA-seq or Spatial-seq (CosMx) data
* Annotation and filtering of cell clusters (scRNA-seq)
* Regional clustering and annotation of cells (Spatial-seq)
* Differential gene expression analysis across regions and/or age groups
* Definition of a set of differentially expressed genes that shift during aging and under heterochronic reconstituiton
* Perform differential cell-cell interaction analysis
* Co-integration Spatial- and scRNA-seq data
* Aging cell neighborhood analysis



## Associated Manuscript
TODO: link to manuscript

# Installation Instructions
Please install the R packages listed in the following 'Dependencies' section.


## Dependencies
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
