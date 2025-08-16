#### CosMx integration + neighborhood (region) analysis ########################
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(future)
  library(dbscan)   # frNN
  library(ggplot2)
})

# ---------- Config ----------
plan(multisession, workers = 20)
options(future.globals.maxSize = 1e10)  # ~10 GB

main_folder   <- "/PATH/TO/MAIN/FOLDER"
input_folder  <- file.path(main_folder, "input_data")
output_folder <- file.path(main_folder, "output_data")
projectID     <- "ReC_wIFNpanel_spatial"

dir.create(input_folder,  recursive = TRUE, showWarnings = FALSE)
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# ---------- Helpers ----------
clear_reductions <- function(x) { x@reductions <- list(); x }

make_numap <- function(obj, x_col = "y_slide_mm", y_col = "x_slide_mm") {
  em <- as.matrix(obj@meta.data[, c(x_col, y_col)])
  colnames(em) <- c("UMAP_1", "UMAP_2")
  obj@reductions$numap <- CreateDimReducObject(embeddings = em, key = "UMAP_", assay = DefaultAssay(obj))
  obj
}

compute_neighborhood_stats <- function(coords_um, labels, radius_um = 100) {
  coords_um <- as.matrix(coords_um)
  lab <- factor(labels)
  nn   <- dbscan::frNN(coords_um, eps = radius_um)
  L    <- levels(lab)
  out  <- matrix(0L, nrow(coords_um), length(L), dimnames = list(NULL, L))
  
  for (i in seq_along(nn$id)) {
    ids <- nn$id[[i]]
    if (length(ids)) {
      keep <- nn$dist[[i]] > 0  # drop self (distance==0)
      if (any(keep)) {
        tab <- table(lab[ids[keep]])
        out[i, names(tab)] <- as.integer(tab)
      }
    }
  }
  # row-wise z-score
  zfun <- function(v) { s <- sd(v); if (s == 0) rep(0, length(v)) else (v - mean(v))/s }
  out  <- t(apply(out, 1, zfun))
  out[is.na(out) | is.infinite(out)] <- 0
  out
}

# ---------- Load CosMx Seurat .rds (one per section) ----------
files <- c(
  Y01 = "/PATH/TO/DATA/Young_01.RDS",
  Y02 = "/PATH/TO/DATA/Young_02.RDS",
  O01 = "/PATH/TO/DATA/Old_01.RDS",
  O02 = "/PATH/TO/DATA/Old_02.RDS"
)

sample_names <- c(
  Y01 = "Young_01",
  Y02 = "Young_02",
  O01 = "Old_01",
  O02 = "Old_02"
)

Spatial_Seq_list <- lapply(names(files), function(k) {
  obj <- readRDS(files[[k]])
  obj$sample <- sample_names[[k]]
  obj <- clear_reductions(obj)
  RenameCells(obj, add.cell.id = k)
})
names(Spatial_Seq_list) <- sample_names

# --- optional manual flips to align sections (edit/remove as needed) ---
# preview (quiet, base ggplot)
lapply(Spatial_Seq_list, \(o) print(ggplot(o@meta.data, aes(y_slide_mm, x_slide_mm))+
   geom_point(size=0.3, alpha=0.6)+theme_void()))

Spatial_Seq_list[["Young_01"]]@meta.data$y_slide_mm <- -Spatial_Seq_list[["Young_01"]]@meta.data$y_slide_mm
Spatial_Seq_list[["Young_02"]]@meta.data$y_slide_mm <- -Spatial_Seq_list[["Young_02"]]@meta.data$y_slide_mm
Spatial_Seq_list[["Old_01"]]  @meta.data$y_slide_mm <- -Spatial_Seq_list[["Old_01"]]  @meta.data$y_slide_mm
Spatial_Seq_list[["Old_02"]]  @meta.data$y_slide_mm <- -Spatial_Seq_list[["Old_02"]]  @meta.data$y_slide_mm

#Confirm orientation
lapply(Spatial_Seq_list, \(o) print(ggplot(o@meta.data, aes(y_slide_mm, x_slide_mm))+
                                      geom_point(size=0.3, alpha=0.6)+theme_void()))

# ----------------------------------------------------------------------

# ---------- Merge & basic cleaning ----------
brain <- Reduce(merge, Spatial_Seq_list)
rm(Spatial_Seq_list); gc()

# drop CosMx-flagged low-quality cells if present
if ("qcCellsFlagged" %in% colnames(brain@meta.data)) {
  brain <- subset(brain, subset = qcCellsFlagged == FALSE)
}

# quick "numap" from slide coordinates for sanity checks
brain <- make_numap(brain, x_col = "y_slide_mm", y_col = "x_slide_mm")
brain$sample <- factor(brain$sample, levels = sample_names, ordered = TRUE)

# ---------- Harmony-only integration ----------
# SCTransform once across all cells (clip keeps dynamic range comparable)


brain[["RNA"]] <- split(brain[["RNA"]], f = brain$sample)

brain <- SCTransform(brain, clip.range = c(-10, 10))
brain <- RunPCA(brain, npcs = 30, verbose = F)

brain <- IntegrateLayers(
  object = brain, normalization.method = "SCT",method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = T,  k.weight = 100
)

brain <- FindNeighbors(brain, reduction = "harmony", dims = 1:30)
brain <- FindClusters(brain, resolution = 0.4, cluster.name = "harmony_clusters")
brain <- RunUMAP(brain, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

(p3 <- DimPlot(
  brain,
  reduction = "umap.harmony",
  group.by = c("harmony_clusters"), split.by = 'sample',
  combine = T
))

DefaultAssay(brain) <- 'RNA'
brain <- JoinLayers(brain)

DefaultAssay(brain) <- 'SCT'
brain <- SCTransform(brain, assay = "RNA", verbose = T, clip.range = c(-10, 10))

# ---------- Neighborhood (region) analysis ----------
# Convert slide mm → µm and compute local composition per cell
brain$x_um <- brain$x_slide_mm * 1000
brain$y_um <- brain$y_slide_mm * 1000
brain$cell_grouping <- paste0("Cluster_", as.character(brain$harmony_clusters))

meta <- brain@meta.data
samples <- levels(brain$sample)
celltypes <- unique(meta$cell_grouping)

nbor_mat <- matrix(0, nrow(meta), length(celltypes))
colnames(nbor_mat) <- celltypes

for (s in samples) {
  idx <- which(meta$sample == s)
  sub_coords <- meta[idx, c("x_um", "y_um")]
  sub_labels <- meta[idx, "cell_grouping", drop = TRUE]
  sub_stats  <- compute_neighborhood_stats(sub_coords, sub_labels, radius_um = 100)
  nbor_mat[idx, colnames(sub_stats)] <- sub_stats
}

# PCA on neighborhood composition (already row-zscored above)
pca_nb <- prcomp(nbor_mat, center = FALSE, scale. = FALSE)
emb_nb <- pca_nb$x

# K-means on PCs → regional clusters (adjust k as needed)
set.seed(345)
k <- 25
km <- kmeans(emb_nb, centers = k, nstart = 20)
brain$regionCluster_kmeans <- factor(km$cluster)

# optional quick looks (no custom themes)
# DimPlot(brain, reduction="umap.harmony", group.by="harmony_clusters", split.by="sample", label=TRUE)
# DimPlot(brain, reduction="numap", group.by="regionCluster_kmeans", split.by="sample")

# ---------- Save ----------
saveRDS(brain, file = file.path(output_folder, paste0(projectID, "_SCT_Harmony_neighborhood.rds")))




