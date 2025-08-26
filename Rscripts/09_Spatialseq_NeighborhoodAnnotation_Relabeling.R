#### CosMx: curate crude region clusters into whole regions ####################
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table) # rbindlist
  library(ggplot2)
})

# --- Config / I/O -------------------------------------------------------------
main_folder <- "/PATH/TO/MAIN/FOLDER"
output_folder <- file.path(main_folder, "output_data")
projectID <- "ReC_wIFNpanel_spatial"

brain <- readRDS(file.path(output_folder, paste0(projectID, "_SCT_Harmony_neighborhood.rds")))

# helpful: ensure a stable per-cell key
if (!"cell_index" %in% colnames(brain@meta.data)) {
  brain$cell_index <- rownames(brain@meta.data)
}

# crude region field from k-means step (rename if yours differs)
# - assumed present as `region_crude`
stopifnot("region_crude" %in% colnames(brain@meta.data))

# a clean color palette for final plots
target_levels <- c("cor", "cer", "cc", "hip", "mid", "plx", "str", "thal", "wm", "vent", "other")
color_code_target_region <- setNames(
  Seurat::DiscretePalette(length(target_levels), palette = "polychrome"),
  target_levels
)

# --- helper: interactive selection then set a label ---------------------------
assign_by_selector <- function(obj, subset_expr, set_to, group_col = "target_region") {
  p <- DimPlot(
    subset(obj, subset = {{ subset_expr }}),
    reduction = "numap", group.by = group_col, label = FALSE, cols = "polychrome"
  )
  picked <- CellSelector(plot = p)
  if (length(picked) > 0) {
    obj@meta.data[picked, group_col] <- set_to
  }
  obj
}

# container to collect per-sample annotations if you iterate multiple samples
meta_data_list <- list()

# ============================ annotate ONE sample =============================
sample_to_analyze <- "Old_02" # <-- change to your sample label
current_brain <- subset(brain, subset = sample == sample_to_analyze)

# 1) initialize a working column from the crude clusters
current_brain$target_region <- "other"
seed_map <- c(
  Thalamus              = "thal",
  Cerebellum            = "cer",
  Cortex                = "cor",
  Striatum              = "str",
  White_matter          = "White_matter", # refined further below
  Hippocampus           = "hip",
  Midbrain_Pons_Medulla = "Midbrain_Pons_Medulla", # refined below
  ChoroidPlexus         = "plx",
  Ventricle             = "vent"
)
idx <- match(current_brain$region_crude, names(seed_map))
keep <- !is.na(idx)
current_brain$target_region[keep] <- unname(seed_map[idx[keep]])

# sanity check view
DimPlot(current_brain, reduction = "numap", group.by = "target_region", cols = "polychrome")

# 2) fold cerebellar **white matter** back into cerebellum
current_brain <- assign_by_selector(
  current_brain,
  subset_expr = region_crude != "Cerebellum",
  set_to      = "cer"
)

# keep CP intact (override if needed)
current_brain$target_region[current_brain$region_crude == "ChoroidPlexus"] <- "plx"

# 3) select **corpus callosum** out of white matter
current_brain <- assign_by_selector(
  current_brain,
  subset_expr = target_region == "White_matter",
  set_to      = "cc"
)

# 4) refine **cortex** within cortex-labeled crude cluster
current_brain <- assign_by_selector(
  current_brain,
  subset_expr = region_crude == "Cortex",
  set_to      = "cor"
)

# 5) refine **hippocampus** from all still-ambiguous areas
current_brain <- assign_by_selector(
  current_brain,
  subset_expr = target_region %in% c("other", "White_matter", "Midbrain_Pons_Medulla"),
  set_to      = "hip"
)
# keep thalamus intact (override if needed)
current_brain$target_region[current_brain$region_crude == "Thalamus"] <- "thal"

# 6) refine **striatum** within striatum-labeled crude cluster
current_brain <- assign_by_selector(
  current_brain,
  subset_expr = region_crude == "Striatum",
  set_to      = "str"
)

# 7) refine **midbrain/pons/medulla** (use target_region subset to avoid overriding new labels)
current_brain <- assign_by_selector(
  current_brain,
  subset_expr = target_region == "Midbrain_Pons_Medulla",
  set_to      = "mid"
)

# 8) finalize white matter
current_brain$target_region[current_brain$target_region == "White_matter"] <- "wm"

# 9) clean up everything else
current_brain$target_region[!(current_brain$target_region %in% target_levels)] <- "other"
current_brain$target_region <- factor(current_brain$target_region, levels = target_levels, ordered = TRUE)

# visualize result for this sample
DimPlot(current_brain,
  reduction = "numap", group.by = "target_region",
  cols = color_code_target_region, raster = FALSE
)

# stash annotated metadata
meta_data_list[[sample_to_analyze]] <- current_brain@meta.data

# ===================== (OPTIONAL) annotate more samples =======================
# repeat the block above for each sample; then bind and map back:

# --- assemble and map back to full object ------------------------------------
meta_data_assembled <- data.table::rbindlist(meta_data_list, fill = TRUE)
meta_data_assembled <- as.data.frame(meta_data_assembled)

map_vec <- setNames(meta_data_assembled$target_region, meta_data_assembled$cell_index)
brain$target_region <- unname(map_vec[brain$cell_index])
brain$target_region[is.na(brain$target_region)] <- "other"
brain$target_region <- factor(brain$target_region, levels = target_levels, ordered = TRUE)

# compare crude vs curated, per sample
# (replace with your preferred palettes)
DimPlot(brain, reduction = "numap", group.by = "region_crude", split.by = "sample")
DimPlot(brain,
  reduction = "numap", group.by = "target_region",
  split.by = "sample", cols = color_code_target_region
)

# save
saveRDS(brain, file = file.path(
  output_folder,
  paste0(projectID, "_Harmony_neighborhood_withRegionLabels.rds")
))
