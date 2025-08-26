#### CosMx immune extraction + ReC co-integration + label transfer #############
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(harmony)
})

set.seed(42)

# -------------------- CONFIG --------------------
proj_dir <- "/PATH/TO/PROJECT"
in_spatial_rds <- file.path(
  proj_dir, "input_data",
  "_Harmony_neighborhood_withRegionLabels_mapmycells_annotated.rds"
)
in_rec_rds <- file.path(proj_dir, "input_data", "ReC_scRNA_annotated.rds") # <-- your ReC object
out_dir <- file.path(proj_dir, "outputs_spatialReC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# fields used later
sample_field <- "sample" # adjust if your per-section tag differs
mapmy_field <- "cellMap_class" # MapMyCells class (e.g., "Immune")
mapmy_prob <- "cellMap_classProb" # MapMyCells class probability
harmony_red <- "umap.harmony" # Harmony UMAP name in your object


# -------------------- UTILITIES --------------------
clear_reductions <- function(seurat_obj) {
  seurat_obj@reductions <- list()
  seurat_obj
}

has_genes <- function(obj, genes) {
  present <- intersect(rownames(obj), genes)
  if (length(present) == 0) {
    warning(
      "None of ", paste(genes, collapse = ", "),
      " are present in this object."
    )
  }
  present
}

integrate_layers_harmony <- function(
    obj,
    batch_key = "sample",
    npcs = 30,
    harmony_dims = 1:30,
    harmony_resolution = 0.4,
    k.weight.harmony = 100,
    clip.range = c(-10, 10),
    verbose = TRUE) {
  stopifnot(batch_key %in% colnames(obj@meta.data))

  # 0) Quick batch summary
  batch_vec <- obj[[batch_key]][, 1]
  nbatches <- length(unique(batch_vec))
  if (verbose) {
    message(sprintf(
      "[integrate_layers_harmony] %s batches in '%s': %s",
      nbatches, batch_key, paste(unique(batch_vec), collapse = ", ")
    ))
  }

  # 1) Split RNA into layers per batch (Seurat v5 layer workflow)
  obj[["RNA"]] <- split(obj[["RNA"]], f = batch_vec)

  # Drop truly empty layers, if any (can happen after filtering)
  nonempty <- vapply(Layers(obj[["RNA"]]), function(ln) ncol(GetAssayData(obj[["RNA"]][[ln]], slot = "counts")) > 0, logical(1))
  if (!all(nonempty)) {
    if (verbose) {
      message(
        "[integrate_layers_harmony] Dropping empty layers: ",
        paste(Layers(obj[["RNA"]])[!nonempty], collapse = ", ")
      )
    }
    obj[["RNA"]] <- obj[["RNA"]][Layers(obj[["RNA"]])[nonempty]]
    # Recompute batch count after dropping empties
    nbatches <- length(Layers(obj[["RNA"]]))
  }

  # 2) Per-layer normalization + PCA
  obj <- SCTransform(obj, clip.range = clip.range, verbose = verbose)
  obj <- RunPCA(obj, npcs = npcs, verbose = verbose)

  # Keep dims within available PCs
  max_pc <- ncol(Embeddings(obj, "pca"))
  dims_umap <- harmony_dims[harmony_dims <= max_pc]
  if (length(dims_umap) == 0) stop("No valid PCs available for downstream steps.")

  # 3) Harmony layer integration (only if >=2 batches)
  if (nbatches >= 2) {
    obj <- IntegrateLayers(
      object = obj,
      normalization.method = "SCT",
      method = HarmonyIntegration,
      orig.reduction = "pca",
      new.reduction = "harmony",
      k.weight = k.weight.harmony,
      verbose = verbose
    )
    red_to_use <- "harmony"
  } else {
    warning(
      "[integrate_layers_harmony] Only one batch/layer detected after split; ",
      "skipping Harmony and proceeding with PCA space."
    )
    red_to_use <- "pca"
  }

  # 4) Graph / cluster / UMAP
  obj <- FindNeighbors(obj, reduction = red_to_use, dims = dims_umap)
  obj <- FindClusters(obj, resolution = harmony_resolution, cluster.name = paste0(red_to_use, "_clusters"))
  obj <- RunUMAP(obj, reduction = red_to_use, dims = dims_umap, reduction.name = paste0("umap.", red_to_use))

  obj
}


# -------------------- LOAD SPATIAL (CosMx) ------------------------------------
spatial_all <- readRDS(file.path(out_dir, "_Harmony_neighborhood_withRegionLabels_mapmycells_annotated.rds"))
DefaultAssay(spatial_all) <- if ("SCT" %in% Assays(spatial_all)) "SCT" else "RNA"

# quick look: immune marker expression across Harmony UMAP
gImm <- c("Tmem119", "Ptprc")
gImm <- has_genes(spatial_all, gImm)
if (length(gImm)) {
  print(FeaturePlot(spatial_all, reduction = harmony_red, features = gImm, order = TRUE))
}

# -------------------- IMMUNE SUBSETTING ---------------------------------------
# Find immune cell cluster - here it's cluster number 13. This will be different for what you're doing!

DefaultAssay(spatial_all) <- "RNA"
immune_ids <- row.names(spatial_all@meta.data %>% filter(harmony_clusters %in% c(13)))

spatial_immune <- subset(spatial_all, cells = immune_ids)
spatial_immune$project <- "CosMx"


# -------------------- RE-INTEGRATE IMMUNE (CosMx only) ------------------------
# (helps tighten immune structure across sections before mixing with scRNA)
DefaultAssay(spatial_immune) <- "RNA"

spatial_immune <- integrate_layers_rpca_harmony(
  obj                = spatial_immune,
  batch_key          = "sample",
  npcs               = 30,
  rpca_dims          = 1:30, # ignored (kept for API compatibility)
  harmony_dims       = 1:30,
  rpca_resolution    = 0.4, # ignored
  harmony_resolution = 0.4,
  k.weight.rpca      = 100, # ignored
  k.weight.harmony   = 100
)

# sanity plots
print(DimPlot(spatial_immune, reduction = "umap.harmony", group.by = "harmony_clusters", label = TRUE))
print(DimPlot(spatial_immune, reduction = "umap.harmony", group.by = sample_field))

saveRDS(spatial_immune, file.path(out_dir, "cosmx_immune_SCT_harmony.rds"))

# -------------------- LOAD ReC scRNA (REFERENCE) ------------------------------
rec <- readRDS(in_rec_rds)
cosmx_im <- spatial_immune
# --- 1) Make them comparable: keep RNA only, add a 'project' tag, align features
DefaultAssay(rec_sc) <- "RNA"
DefaultAssay(cosmx_im) <- "RNA"

rec_sc <- NormalizeData(rec_sc)
cosmx_im <- NormalizeData(cosmx_im)

rec_sc <- DietSeurat(rec_sc, assays = "RNA")
rec_sc$project <- "scRNA"
cosmx_im <- DietSeurat(cosmx_im, assays = "RNA")
cosmx_im$project <- "CosMx"

common_features <- intersect(rownames(rec_sc), rownames(cosmx_im))
rec_sc <- subset(rec_sc, features = common_features)
cosmx_im <- subset(cosmx_im, features = common_features)

# --- 2) Merge -> immune_mix (Seurat v5 layers come next inside the helper)
immune_mix <- merge(rec_sc, y = cosmx_im)
immune_mix <- JoinLayers(immune_mix) # important for Seurat v5 layer ops

# --- 3) Layer-based Harmony across 'project' (scRNA vs CosMx) ---------------
# (Assumes you have integrate_layers_harmony() defined from the previous step.)
immune_mix <- integrate_layers_harmony(
  obj                = immune_mix,
  batch_key          = "project", # << key point
  npcs               = 30,
  harmony_dims       = 1:30,
  harmony_resolution = 0.4,
  k.weight.harmony   = 100,
  clip.range         = c(-10, 10),
  verbose            = TRUE
)

# --- 4) Quick sanity checks ---------------------------------------------------
# Co-embedding by project and cluster
DimPlot(immune_mix,
  reduction = "umap.harmony",
  group.by = c("project", "harmony_clusters")
)

# Microglia / immune markers across the co-embedding
FeaturePlot(immune_mix,
  reduction = "umap.harmony",
  features = c("Tmem119", "P2ry12", "Ptprc", "Spp1", "Gpnmb"),
  order = TRUE
)


# --- 5) Manual label transfer scaffold (no auto TransferData) ----------------
# Visualize known scRNA labels (e.g., rec_sc$cell_type) and use them to guide manual calls
DimPlot(subset(immune_mix, subset = project == "scRNA"),
  reduction = "umap.harmony", group.by = "cell_type", label = TRUE
)

# Create your manual mapping from CosMx Harmony clusters -> chosen microglia states
# (Fill these in by inspection)
manual_map <- c(
  `0` = "Homeostatic",
  `1` = "Interferon_responsive",
  `2` = "MacrophagePeripheral",
  `3` = "MyelinPositive"
  # ...add the rest as you decide
)

cosmx_only <- subset(immune_mix, subset = project == "CosMx")
cosmx_only$manual_label <- manual_map[as.character(cosmx_only$harmony_clusters)]

# Put manual labels back onto the full object
immune_mix$manual_label <- NA_character_
immune_mix$manual_label[colnames(cosmx_only)] <- cosmx_only$manual_label

# Visualize your manual labels on the co-embedding
DimPlot(immune_mix, reduction = "umap.harmony", group.by = "manual_label")




## --------------------------------------------
## Map back scRNA-seq–guided labels to spatial
## --------------------------------------------

stopifnot(inherits(spatial_all, "Seurat"))
stopifnot(inherits(immune_mix, "Seurat"))
stopifnot("project" %in% colnames(immune_mix@meta.data))

# 1) Decide which label to use as the "cell type" coming from immune_mix
label_type_col <- dplyr::case_when(
  "manual_label" %in% colnames(immune_mix@meta.data) ~ "manual_label", # preferred: your manual transfer
  "ReC_label" %in% colnames(immune_mix@meta.data) ~ "ReC_label", # fallback: TransferData results
  TRUE ~ "harmony_clusters" # last resort: numeric clusters
)

# 2) If you already created a manual class column, we’ll use it; otherwise derive a coarse class
has_manual_class <- "manual_class" %in% colnames(immune_mix@meta.data)
derive_class <- function(x) {
  # x: vector of type labels (character)
  out <- rep("other", length(x))
  out[grepl("Tcell|T cell|Cd3", x, ignore.case = TRUE)] <- "Tcell"
  out[grepl("NK", x, ignore.case = TRUE)] <- "NKcell"
  out[grepl("Macrophage|Mono", x, ignore.case = TRUE)] <- "MacrophagePeripheral"
  out[grepl("Microglia|Homeostatic|IFN|Interferon|Cytokine|DAM|Myelin",
    x,
    ignore.case = TRUE
  )] <- "Microglia"
  out[grepl("Reconstituted|ReC", x, ignore.case = TRUE)] <- "Reconstituted"
  out
}

# 3) Pull CosMx cells from the co-embedding (immune_mix)
cosmx_ids <- rownames(immune_mix@meta.data)[immune_mix$project == "CosMx"]

# Build a small mapping table (cell_id -> type/class) from immune_mix
type_vec <- immune_mix@meta.data[cosmx_ids, label_type_col, drop = TRUE]
class_vec <- if (has_manual_class) {
  immune_mix@meta.data[cosmx_ids, "manual_class", drop = TRUE]
} else {
  derive_class(as.character(type_vec))
}

map_df <- data.frame(
  cell_id = cosmx_ids,
  map_type = as.character(type_vec),
  map_class = as.character(class_vec),
  stringsAsFactors = FALSE,
  row.names = cosmx_ids
)

# 4) Initialize new columns on the full spatial object as "other"
spatial_all$scRNAseq_cell_type <- "other"
spatial_all$scRNAseq_cell_class <- "other"

# 5) Transfer labels only for overlapping CosMx cells
overlap_ids <- intersect(colnames(spatial_all), rownames(map_df))
spatial_all$scRNAseq_cell_type[overlap_ids] <- map_df[overlap_ids, "map_type"]
spatial_all$scRNAseq_cell_class[overlap_ids] <- map_df[overlap_ids, "map_class"]

# 6) Quick sanity plots (use numap if you built it; else Harmony UMAP if present)
red_to_use <- if ("numap" %in% names(spatial_all@reductions)) {
  "numap"
} else if ("umap.harmony" %in% names(spatial_all@reductions)) {
  "umap.harmony"
} else {
  NULL
}

if (!is.null(red_to_use)) {
  print(DimPlot(spatial_all, reduction = red_to_use, group.by = "scRNAseq_cell_type"))
  print(DimPlot(spatial_all, reduction = red_to_use, group.by = "scRNAseq_cell_class"))
}

# Optional: check counts
table(spatial_all$scRNAseq_cell_class)


saveRDS(spatial_all, file.path(out_dir, "_Harmony_neighborhood_withRegionLabels_mapmycells_annotated_ReCcluster_annotated.rds"))
