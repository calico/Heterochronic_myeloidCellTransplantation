#### 0) Setup ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(future)
})

# Reusable helper --------------------------------------------------------------

set.seed(42)
plan(multisession)
options(future.globals.maxSize = 3e9)

# Project structure (set your working directory to the repo root before running)
project_dir <- '/data/novaseq-2024/240523_A01059_0383_BHFNN3DSXC/REQ_24_089_090_cellranger_with_featureCapture/'
input_dir   <- file.path(project_dir, "input_data")         # contains per-sample Cell Ranger folders
output_dir  <- file.path(project_dir, "outputs_scRNA")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#### 1) Discover 10x inputs ----------------------------------------------------
# Expect: input_data/<sample_id>/outs/filtered_feature_bc_matrix/
sample_roots <- list.dirs(input_dir, full.names = TRUE, recursive = FALSE)
matrix_dirs  <- file.path(sample_roots, "outs", "filtered_feature_bc_matrix")
matrix_dirs  <- matrix_dirs[dir.exists(matrix_dirs)]

sample_ids <- basename(sample_roots[dir.exists(file.path(sample_roots, "outs"))])

#### 2) Load samples -----------------------------------------------------------
sc_list <- vector("list", length(matrix_dirs))
names(sc_list) <- sample_ids

for (i in seq_along(matrix_dirs)) {
  fbm <- Read10X(data.dir = matrix_dirs[i])
  # Create RNA assay from "Gene Expression"
  obj <- CreateSeuratObject(
    counts      = fbm[["Gene Expression"]],
    min.cells   = 3,
    min.features= 200
  )
  # Add ADT assay from "Antibody Capture" (if present)
  if ("Antibody Capture" %in% names(fbm)) {
    obj@assays$ADT <- CreateAssay5Object(counts = fbm[["Antibody Capture"]])
  }
  obj$sample_id <- names(sc_list)[i]
  sc_list[[i]] <- RenameCells(obj, add.cell.id = obj$sample_id)
  
  
}

#### 3) Merge ------------------------------------------------------------------
sc_object <- Reduce(function(x, y) merge(x, y), sc_list)
rm(sc_list); gc()

#### 4) QC metrics -------------------------------------------------------------
sc_object[["percent_mt"]]  <- PercentageFeatureSet(sc_object, pattern = "^mt-")
sc_object[["percent_rpl"]] <- PercentageFeatureSet(sc_object, pattern = "^rp[sl]")

# Optional EGFP (set to a gene name if applicable, else leave NULL)
sc_object[["percent_egfp"]] <- PercentageFeatureSet(sc_object, pattern = 'NLS-Cas9-NLS-P2A-EGFP')
sc_object$egfp_pos <- ifelse(sc_object$percent_egfp > 0, "Pos", "Neg")

#### 5) (Your) metadata input --------------------------------------------------
# Supply your sample annotations here. Two options:
# (A) Provide a CSV at input_data/metadata.csv with at least: sample_id,tissue,treatment
# (B) Or edit the inline data.frame below.

metadata_path <- file.path(input_dir, "metadata.csv")

if (file.exists(metadata_path)) {
  md <- utils::read.csv(metadata_path, stringsAsFactors = FALSE)
} else {
  md <- data.frame(
    sample_id = unique(sc_object$sample_id),
    tissue    = NA_character_,   # <-- fill in (e.g., "Cortex", "Cerebellum")
    treatment = NA_character_,   # <-- fill in (e.g., "Young", "Old")
    replicate = NA_character_,   # <-- optional
    stringsAsFactors = FALSE
  )
  # write a template so users know the expected format
  utils::write.csv(md, file.path(output_dir, "metadata_template.csv"), row.names = FALSE)
  message("No metadata.csv found. A template was written to outputs_scRNA/metadata_template.csv")
}

# Join metadata by sample_id (no assumptions about naming)
sc_object@meta.data <- sc_object@meta.data %>%
  left_join(md, by = "sample_id")

Idents(sc_object) <- "sample_id"

#### 6) Save raw merged object -------------------------------------------------
saveRDS(sc_object, file.path(output_dir, "sc_object_raw.rds"))

#### 7) QC filtering (generic thresholds) -------------------------------------
qc <- with(sc_object@meta.data,
           nFeature_RNA > 200 & nFeature_RNA < 4000 & percent_mt < 10)
sc_object <- subset(sc_object, cells = colnames(sc_object)[qc])

saveRDS(sc_object, file.path(output_dir, "sc_object_qc_filtered.rds"))

#### 8) Minimal QC plots -------------------------------------------------------
dir.create(file.path(output_dir, "qc_plots"), showWarnings = FALSE)

p_vln <- VlnPlot(sc_object, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"),
                 group.by = "sample_id", pt.size = 0.1) + NoLegend()

#### 9) Normalize, PCA ---------------------------------------------------------
sc_object <- SCTransform(sc_object, verbose = FALSE)
sc_object <- RunPCA(sc_object, npcs = 30, verbose = FALSE)

#### 10) Harmony integration ---------------------------------------------------



sc_object[["RNA"]] <- split(sc_object[["RNA"]], f = sc_object$sample)

sc_object <- SCTransform(sc_object, clip.range = c(-10, 10))
sc_object <- RunPCA(sc_object, npcs = 30, verbose = F)

sc_object <- IntegrateLayers(
  object = sc_object, normalization.method = "SCT",method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = T,  k.weight = 100
)

sc_object <- FindNeighbors(sc_object, reduction = "harmony", dims = 1:30)
sc_object <- FindClusters(sc_object, resolution = 0.4, cluster.name = "harmony_clusters")
sc_object <- RunUMAP(sc_object, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

(p3 <- DimPlot(
  sc_object,
  reduction = "umap.harmony",
  group.by = c("harmony_clusters"), split.by = 'sample',
  combine = T
))

DefaultAssay(sc_object) <- 'RNA'
sc_object <- JoinLayers(sc_object)

#### 11) UMAP panels (generic facets) -----------------------------------------

DimPlot(sc_object, reduction = "umap_harmony", group.by = "tissue")
DimPlot(sc_object, reduction = "umap_harmony", group.by = "treatment")
DimPlot(sc_object, reduction = "umap_harmony", group.by = "sample_id")
DimPlot(sc_object, reduction = "umap_harmony", group.by = "harmony_clusters", label = TRUE)

#### 12) Normalize and scale RNA layer for DEG analysis and plotting -----------------------------------------

#Run preprocessing with conventional RNA normalization and PCA
sc_object <- NormalizeData(sc_object)
sc_object <- FindVariableFeatures(sc_object)
sc_object <- ScaleData(sc_object)
sc_object <- RunPCA(sc_object)

#Only if ADT/Cite-seq data is present
DefaultAssay(sc_object) <- 'ADT'
sc_object <- NormalizeData(sc_object, normalization.method = "CLR", margin = 2)


#### 13) Save processed object -------------------------------------------------
saveRDS(sc_object, file.path(output_dir, "sc_object_harmony.rds"))

#### 14) Session info ----------------------------------------------------------
sessionInfo()