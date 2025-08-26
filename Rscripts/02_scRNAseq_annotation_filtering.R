#### 0) Setup ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(plyr)
})


project_dir <- getwd()
out_dir <- file.path(project_dir, "outputs_scRNA")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

obj_path <- file.path(out_dir, "sc_object_harmony.rds")

sc_object <- readRDS(obj_path)


sc_object$age <- factor(sc_object$age, levels = c("Young", "Old"), ordered = T)


Idents(sc_object) <- "harmony_clusters"

#### 1) Marker discovery per cluster ------------------------------------------
# If you prefer RNA-level markers, set assay_to_use <- "SCT"
DefaultAssay(sc_object) <- "RNA"

markers_by_cluster <- FindAllMarkers(
  object = sc_object,
  only.pos = TRUE,
  min.pct = 0.15,
  logfc.threshold = 0.15,
  assay = "RNA",
  verbose = TRUE
)

saveRDS(markers_by_cluster, file.path(out_dir, "markers_by_cluster.rds"))

#### 2) Manual cluster→cell-type mapping --------------------------------------
# Inspect `markers_by_cluster` and define your mapping below.
# Use your study’s labels: HM, IRM, RM, CRM, Apoe_high, BAM, T/NK, Granulocyte, Ependymal, etc.
# EXAMPLE (edit to your cluster IDs and labels):
cluster_to_celltype <- c(
  `0` = "HM",
  `1` = "IRM",
  `2` = "RM",
  `3` = "CRM",
  `4` = "Apoe_high",
  `5` = "Monocyte",
  `6` = "T/NK",
  `7` = "Granulocyte",
  `8` = "Ependymal",
  `9` = "Microglia"
)

cluster_levels <- levels(sc_object$harmony_clusters)
cluster_to_celltype <- setNames(rep(NA_character_, length(cluster_levels)), cluster_levels) # <- fill in

# apply mapping (will warn if any NA remain)
sc_object$cell_type <- unname(cluster_to_celltype[as.character(sc_object$harmony_clusters)])
if (anyNA(sc_object$cell_type)) {
  warning("Some clusters lack a cell type label. Fill in `cluster_to_celltype` for all cluster IDs.")
}

# We will also establish a less fine-grained cell annotation where we put all Reconstituted cell states as one 'cell class'
sc_object$cell_class <- plyr::mapvalues(
  x    = sc_object$cell_type,
  from = c(rec_like, "Monocyte", "T/NK", "Granulocyte", "Ependymal", "Microglia"),
  to   = c(rep("ReC", length(rec_like)), "Monocyte", "T/NK", "Granulocyte", "Ependymal", "Microglia")
)

# set identities to cell type for plotting
Idents(sc_object) <- "cell_type"

#### 3) Quick figure panels ----------------------------------------------------
# UMAP by cell type and by age
p_ct <- DimPlot(sc_object, reduction = "umap_harmony", group.by = "cell_type", label = TRUE)
p_age <- DimPlot(sc_object, reduction = "umap_harmony", group.by = "age")

# Example gene panels (edit features as useful for your labels)
feature_markers <- c("P2ry12", "Cx3cr1", "Ifit3", "Isg15", "Rplp0", "Rps27a", "Tnf", "Il1b", "Apoe", "Apoe") # replace as needed
feature_markers <- feature_markers[feature_markers %in% rownames(sc_object)]
if (length(feature_markers)) {
  g1 <- FeaturePlot(sc_object, features = feature_markers, reduction = "umap_harmony", combine = TRUE)
}


#### 4) Optional filtering (user-controlled) -----------------------------------
# List any labels to drop (e.g., doublets/unknown/debris).
drop_labels <- c("Doublet", "Unknown", "Debris") # <- edit / keep empty to skip
if (length(drop_labels)) {
  keep <- is.na(sc_object$cell_type) | !(sc_object$cell_type %in% drop_labels)
  sc_object <- subset(sc_object, cells = colnames(sc_object)[keep])
}

#### 5) Save annotated objects -------------------------------------------------
saveRDS(sc_object, file.path(out_dir, "sc_object_annotated.rds"))

# Convenience: a “clean” object with only labeled cell types (no NA)
sc_object_clean <- subset(sc_object, subset = !is.na(cell_type))
saveRDS(sc_object_clean, file.path(out_dir, "sc_object_annotated_clean.rds"))

#### 6) session info -----------------------------------------------------------
sessionInfo()
