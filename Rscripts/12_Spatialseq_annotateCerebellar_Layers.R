suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(RANN) # nn2
  library(data.table)
  library(ggplot2)
})

# --------- Assumptions / columns on spatial_all@meta.data --------------------
# - coordinates in micrometers: cols "x" and "y"
#     (if missing, we’ll create them from x_slide_mm / y_slide_mm * 1000)
# - region field: target_region == "cer" marks cerebellum
# - immune labels mapped back: scRNAseq_cell_class == "Reconstituted" for ReCs
# - MapMyCells annotations (class/subclass) present:
#     cellMap_class: e.g. "28 CB GABA", "29 CB Glut"
#     cellMap_subclass: e.g. "313 CBX Purkinje Gaba", "314 CB Granule Glut"

#-------------------- CONFIG --------------------
proj_dir <- "/PATH/TO/PROJECT"
in_spatial_rds <- file.path(
  proj_dir, "input_data",
  "_Harmony_neighborhood_withRegionLabels_mapmycells_annotated_ReCcluster_annotated.rds"
)
in_rec_rds <- file.path(proj_dir, "input_data", "ReC_scRNA_annotated.rds") # <-- your ReC object
out_dir <- file.path(proj_dir, "outputs_spatialReC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)




# If x/y in µm are missing, create them from slide mm columns
ensure_xy_um <- function(obj, x_col = "x", y_col = "y",
                         x_mm = "y_slide_mm", y_mm = "x_slide_mm") {
  md <- obj@meta.data
  if (!all(c(x_col, y_col) %in% colnames(md))) {
    stopifnot(all(c(x_mm, y_mm) %in% colnames(md)))
    obj@meta.data[[x_col]] <- obj@meta.data[[x_mm]] * 1000
    obj@meta.data[[y_col]] <- obj@meta.data[[y_mm]] * 1000
  }
  obj
}

# Core classifier: given a per-sample data.frame, attach cerebellar layer to ReCs
classify_cerebellar_microglia <- function(df,
                                          dist_um = 30,
                                          type_field = "scRNAseq_cell_class",
                                          rec_value = "Reconstituted",
                                          mm_class_field = "cellMap_class",
                                          mm_subclass_field = "cellMap_subclass",
                                          x_col = "x", y_col = "y") {
  stopifnot(all(c(type_field, mm_class_field, mm_subclass_field, x_col, y_col) %in% colnames(df)))

  # Reconstituted microglia
  df_micro <- df %>% filter(.data[[type_field]] == rec_value)
  if (nrow(df_micro) == 0) {
    return(NULL)
  }

  # Purkinje and granule neuron pools (Allen MapMyCells)
  df_purk <- df %>% filter(
    .data[[mm_class_field]] == "28 CB GABA",
    .data[[mm_subclass_field]] == "313 CBX Purkinje Gaba"
  )
  df_gran <- df %>% filter(
    .data[[mm_class_field]] == "29 CB Glut",
    .data[[mm_subclass_field]] == "314 CB Granule Glut"
  )

  coords_micro <- as.matrix(df_micro[, c(x_col, y_col)])
  dist_to_purk <- rep(Inf, nrow(df_micro))
  dist_to_gran <- rep(Inf, nrow(df_micro))

  if (nrow(df_purk) > 0) {
    nn <- nn2(data = as.matrix(df_purk[, c(x_col, y_col)]), query = coords_micro, k = 1)
    dist_to_purk <- as.numeric(nn$nn.dists[, 1])
  }
  if (nrow(df_gran) > 0) {
    nn <- nn2(data = as.matrix(df_gran[, c(x_col, y_col)]), query = coords_micro, k = 1)
    dist_to_gran <- as.numeric(nn$nn.dists[, 1])
  }

  # Rules:
  # - neither within dist_um -> white_matter
  # - both within -> choose closer (purk -> ML, gran -> GL)
  # - only purk  -> ML
  # - only gran  -> GL
  layer <- pmap_chr(list(dist_to_purk, dist_to_gran), function(dp, dg) {
    in_p <- dp < dist_um
    in_g <- dg < dist_um
    if (!in_p && !in_g) {
      return("white_matter")
    }
    if (in_p && in_g) {
      return(if (dp < dg) "molecular_layer" else "granule_layer")
    }
    if (in_p && !in_g) {
      return("molecular_layer")
    }
    # /* (!in_p &&  in_g) */ 
    return("granule_layer")
  })

  out <- df_micro %>%
    mutate(
      cerebellar_layer = layer,
      dist_to_purkinje = dist_to_purk,
      dist_to_granule = dist_to_gran
    ) %>%
    select(cell_id = any_of(c("cell_id", "cell_index", "barcode", "CellID"))[1], everything())
  if (!"cell_id" %in% colnames(out)) {
    out$cell_id <- rownames(out) # fall back to rownames if no explicit ID column
  }
  out[, c("cell_id", "cerebellar_layer", "dist_to_purkinje", "dist_to_granule")]
}

# ---------- RUN (per sample) and map back to spatial_all ---------------------

# 1) Ensure µm coordinates exist
spatial_all <- ensure_xy_um(spatial_all, x_col = "x", y_col = "y")

# 2) Work only inside the cerebellum region
stopifnot("target_region" %in% colnames(spatial_all@meta.data))
cer_cells <- rownames(spatial_all@meta.data)[spatial_all$target_region == "cer"]
cer_df <- spatial_all@meta.data[cer_cells, , drop = FALSE]
cer_df$cell_id <- rownames(cer_df)

# 3) Split by section if desired (recommended to keep distances local)
by_sample <- split(cer_df, cer_df$sample)

res_list <- lapply(by_sample, function(dd) {
  classify_cerebellar_microglia(
    df = dd,
    dist_um = 30, # tweak if your CosMx pixel->µm differs
    type_field = "scRNAseq_cell_class", # ReC label you mapped back
    rec_value = "Reconstituted",
    mm_class_field = "cellMap_class",
    mm_subclass_field = "cellMap_subclass",
    x_col = "x", y_col = "y"
  )
})

# Bind results (some samples may return NULL if no ReCs inside 'cer')
res_bind <- data.table::rbindlist(Filter(Negate(is.null), res_list))
res_bind <- as.data.frame(res_bind)
rownames(res_bind) <- res_bind$cell_id

# 4) Write back to spatial_all metadata (initialize others as NA)
spatial_all$cerebellar_layer <- NA_character_
overlap <- intersect(colnames(spatial_all), rownames(res_bind))
spatial_all$cerebellar_layer[overlap] <- res_bind[overlap, "cerebellar_layer"]

# 5) Quick counts
print(table(spatial_all$cerebellar_layer, useNA = "ifany"))

# 6) Quick spatial visualization (uses numap if present; otherwise Harmony UMAP)
red_to_use <- if ("numap" %in% names(spatial_all@reductions)) "numap" else if ("umap.harmony" %in% names(spatial_all@reductions)) "umap.harmony" else NULL
if (!is.null(red_to_use)) {
  print(DimPlot(spatial_all, reduction = red_to_use, group.by = "cerebellar_layer"))
}

# 7) Optional: raw XY plot for one sample (overlay Purkinje / granule landmarks)
one <- names(by_sample)[1]
plot_df <- by_sample[[one]]
plot_layers <- res_list[[one]]
if (!is.null(plot_layers)) {
  gg <- ggplot() +
    geom_point(
      data = plot_df %>% dplyr::filter(
        cellMap_class == "28 CB GABA",
        cellMap_subclass == "313 CBX Purkinje Gaba"
      ),
      aes(x = x, y = y), color = "hotpink", alpha = 0.25, size = 0.4
    ) +
    geom_point(
      data = plot_df %>% dplyr::filter(
        cellMap_class == "29 CB Glut",
        cellMap_subclass == "314 CB Granule Glut"
      ),
      aes(x = x, y = y), color = "purple", alpha = 0.15, size = 0.4
    ) +
    geom_point(data = plot_layers, aes(x = x, y = y, color = cerebellar_layer), size = 0.5) +
    theme_void() +
    coord_equal()
  print(gg)
}



saveRDS(spatial_all, file.path(out_dir, "_Harmony_neighborhood_withRegionLabels_mapmycells_annotated_ReCcluster_annotated_CerebellarLayer_annotated.rds"))
