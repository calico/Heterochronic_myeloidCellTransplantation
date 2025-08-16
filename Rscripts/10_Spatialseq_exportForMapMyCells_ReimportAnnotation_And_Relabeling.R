#### MapMyCells export/import for CosMx (or other) Seurat objects ##############
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  # choose one of the two export routes below:
  # 1) anndata (recommended)  2) SeuratDisk
  library(anndata)
  # library(SeuratDisk)
})

set.seed(1)

# -------------------- CONFIG --------------------
main_dir     <- "/PATH/TO/PROJECT"
input_rds    <- file.path(main_dir, "input_data", "_Harmony_neighborhood_withRegionLabels.rds")
out_dir      <- file.path(main_dir, "mapmycells")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# which assay to export to MapMyCells (SCT or RNA)
export_assay <- "SCT"     # change to "RNA" if desired

# path to write the .h5ad that you upload to MapMyCells
out_h5ad     <- file.path(out_dir, "spatial_integrated_counts.h5ad")

# -------------------- LOAD ----------------------
obj <- readRDS(input_rds)
DefaultAssay(obj) <- export_assay

# -------------------- EXPORT to .h5ad (upload to MapMyCells) ------------------
# MapMyCells expects cells x genes counts (sparse is fine).
# We'll export the assay data as-is (transpose to cells x genes).
message(sprintf("Exporting '%s' assay to h5ad: %s", export_assay, out_h5ad))
mat <- t(GetAssayData(obj, assay = export_assay, slot = "data"))
# Prefer counts if available:
if ("counts" %in% Assays(obj)[[export_assay]]@layers) {
  mat <- t(GetAssayData(obj, assay = export_assay, slot = "counts"))
}

# coerce to sparse if needed
mat <- as(mat, "dgCMatrix")

# ensure gene names unique
genes <- make.unique(colnames(mat))
colnames(mat) <- genes
cells <- rownames(mat)

ad <- AnnData(
  X   = mat,
  var = data.frame(gene_symbol = genes, row.names = genes, check.names = FALSE),
  obs = data.frame(cell_id = cells, row.names = cells, check.names = FALSE)
)
write_h5ad(ad, out_h5ad)
message("Done. Upload the .h5ad to MapMyCells and download the CSV report when ready.")

# -------------------- IMPORT MapMyCells report & annotate ---------------------
# The MMC CSV has a 4-line header before the table in most releases; if your
# file differs, set `skip_lines` accordingly.
# AFTER you run MapMyCells online, download the CSV report and point here:
mmc_csv      <- file.path(out_dir, "spatial_integrated_counts_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping.csv")

skip_lines <- 4
stopifnot(file.exists(mmc_csv))

mmc <- suppressMessages(read_csv(mmc_csv, skip = skip_lines))
mmc <- as.data.frame(mmc)

required_cols <- c("cell_id","class_name","class_bootstrapping_probability",
                   "subclass_name","subclass_bootstrapping_probability")
missing_cols  <- setdiff(required_cols, colnames(mmc))
if (length(missing_cols)) {
  stop("MapMyCells CSV is missing columns: ", paste(missing_cols, collapse = ", "))
}

# Harmonize IDs and join onto your Seurat metadata
meta <- obj@meta.data
meta$cell_id <- rownames(meta)

annot <- mmc |>
  transmute(
    cell_id,
    cellMap_class        = class_name,
    cellMap_classProb    = as.numeric(class_bootstrapping_probability),
    cellMap_subclass     = subclass_name,
    cellMap_subclassProb = as.numeric(subclass_bootstrapping_probability)
  )

meta2 <- meta |> left_join(annot, by = "cell_id")
rownames(meta2) <- meta2$cell_id
obj@meta.data <- meta2[, setdiff(colnames(meta2), "cell_id")]

# quick match sanity check
match_rate <- mean(!is.na(obj$cellMap_class))
message(sprintf("MapMyCells annotations matched for %.1f%% of cells.",
                100 * match_rate))

# -------------------- QUICK LOOKS --------------------
p1 <- DimPlot(obj, reduction = "umap.harmony", group.by = "cellMap_class", label = FALSE) +
  ggtitle("MapMyCells class")
p2 <- DimPlot(obj, reduction = "umap.harmony", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("Harmony clusters")
p1 + p2

# confusion / co-occurrence (clusters × class)
tab <- table(obj$seurat_clusters, obj$cellMap_class)
print(pheatmap::pheatmap(tab, cluster_rows = TRUE, cluster_cols = TRUE))

# save annotated object
saveRDS(obj, file.path(out_dir, "_Harmony_neighborhood_withRegionLabels_mapmycells_annotated.rds"))
