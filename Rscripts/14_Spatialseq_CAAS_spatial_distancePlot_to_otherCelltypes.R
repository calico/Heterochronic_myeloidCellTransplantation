#### Spatial proximity of signature scores vs. nearby cell types ################
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(RANN) # fast nearest neighbors
  library(parallel)
  library(ggplot2)
})

set.seed(1)

#-------------------- CONFIG --------------------
proj_dir <- "/PATH/TO/PROJECT"
in_rds <- file.path(
  "/PATH/TO/output",
  "_Harmony_neighborhood_withRegionLabels_mapmycells_annotated_ReCcluster_annotated_CerebellarLayer_annotated_scScore_quantified.rds"
)
out_dir <- "/PATH/TO/output"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# Coordinates in metadata (mm). Set convert_mm_to_um=TRUE if you want µm.
x_col <- "y_slide_mm"
y_col <- "x_slide_mm"
convert_mm_to_um <- TRUE

# Which cells are your “microglia/immune” of interest?
immune_class_col <- "immune_class"
immune_class_val <- "Reconstituted" # e.g., “Reconstituted”

# Which per-cell signature to test? Make sure you have run the score calculation script on this Spatial-seq object first
signature_col <- "Aging_Cerebellum_shared_WT_ReC"

# Which metadata fields hold sample→age/replicate? Provide either mappings or
# leave NULL if you already have clean columns in the object.
sample_col <- "sample"
age_map <- c(Young_01 = "young", Young_02 = "young", Old_01 = "old", Old_02 = "old")
rep_map <- c(Young_01 = "rep1", Young_02 = "rep2", Old_01 = "rep1", Old_02 = "rep2")

# Region filter (optional): set to NULL to use all cells, or e.g. "cer"
region_col <- "target_region"
region_keep <- NULL # e.g., "cer" to restrict to cerebellum

# Which label column defines “cell types” for neighbors?
type_col <- "cellMap_subclass" # often “cellMap_subclass” (MapMyCells subclass)
# For immune cells, optionally override that label with your immune annotations:
immune_override_col <- "immune_class" # set to NULL to disable override

# Distance–binning parameters (units follow your coords after conversion)
bin_max <- 80 # µm
bin_step <- 1 # µm
bin_window <- 30 # µm  (sliding window width)
cores <- max(1, parallel::detectCores() - 1)

# ---------- Load & prepare ----------
spatial_all <- readRDS(in_rds)

md <- spatial_all@meta.data %>%
  mutate(cell_id = rownames(.))


# Coordinates
stopifnot(x_col %in% colnames(md), y_col %in% colnames(md))
md <- md %>%
  transmute(
    cell_id,
    x = .data[[x_col]],
    y = .data[[y_col]],
    across(-c(cell_id), identity, .names = "{.col}") # keep other columns
  )

if (convert_mm_to_um) {
  md$x <- md$x * 1000
  md$y <- md$y * 1000
}

# Build a “celltype” field
stopifnot(type_col %in% colnames(md))
md$celltype <- md[[type_col]]

# Optionally override celltype for immune cells so they appear as “Microglia/Tcell/MΦ…”
if (!is.null(immune_override_col) && immune_override_col %in% colnames(md)) {
  idx <- !is.na(md[[immune_override_col]]) & md[[immune_override_col]] != "other"
  md$celltype[idx] <- md[[immune_override_col]][idx]
}

# Age/replicate
if (!is.null(age_map)) {
  stopifnot(sample_col %in% colnames(md))
  samp <- sub("^.*_", "", md[[sample_col]])
  # robust map by detecting keys in sample string
  pick <- function(x, map) {
    hit <- rep(NA_character_, length(x))
    for (k in names(map)) {
      hit[grepl(k, md[[sample_col]], fixed = TRUE)] <- map[[k]]
    }
    hit
  }
  md$age <- pick(md[[sample_col]], age_map)
  md$replicate <- pick(md[[sample_col]], rep_map)
  md$age[is.na(md$age)] <- if ("treatment" %in% colnames(md)) md$treatment[is.na(md$age)] else "unknown"
  md$replicate[is.na(md$replicate)] <- if ("replicate" %in% colnames(md)) md$replicate[is.na(md$replicate)] else "rep"
} else {
  if (!("age" %in% colnames(md))) md$age <- "all"
  if (!("replicate" %in% colnames(md))) md$replicate <- "rep"
}

# Choose immune/microglia of interest & signature
stopifnot(signature_col %in% colnames(md))
stopifnot(immune_class_col %in% colnames(md))

mg <- md %>% filter(.data[[immune_class_col]] == immune_class_val)
if (nrow(mg) == 0) stop("No immune cells found that match: ", immune_class_val)

# Normalize signature to [0,1] separately per (age,replicate) to reduce batch effects
mg <- mg %>%
  group_by(age, replicate) %>%
  mutate(
    sig_raw = .data[[signature_col]],
    sig_norm = {
      rng <- range(sig_raw, na.rm = TRUE)
      if (diff(rng) == 0) rep(0, length(sig_raw)) else (sig_raw - rng[1]) / diff(rng)
    }
  ) %>%
  ungroup()

# Define target cell types (neighbors)
# keep reasonably frequent types and add canonical immune types explicitly
counts <- sort(table(md$celltype), decreasing = TRUE)
frequent <- names(counts[counts >= 200])
target_types <- setdiff(unique(c(frequent, "ReC")), immune_class_val)
target_types <- intersect(target_types, unique(md$celltype))
if (length(target_types) == 0) stop("No target cell types found; check `type_col` and thresholds.")

message("Target types: ", paste(target_types, collapse = ", "))

# ---------- Distance calc helpers ----------
kd_min_dist <- function(query_xy, ref_xy) {
  # nn2 returns Euclidean distances; k=1 is min distance
  out <- nn2(data = ref_xy, query = query_xy, k = 1)$nn.dists[, 1]
  as.numeric(out)
}

compute_min_dists_for_group <- function(mg_grp, all_cells, target_types) {
  query <- as.matrix(mg_grp[, c("x", "y")])
  res <- matrix(NA_real_, nrow(mg_grp), length(target_types))
  colnames(res) <- paste0("dist_to_", target_types)

  for (j in seq_along(target_types)) {
    ref <- all_cells %>% filter(celltype == target_types[j])
    if (nrow(ref) > 0) {
      ref_xy <- as.matrix(ref[, c("x", "y")])
      res[, j] <- kd_min_dist(query, ref_xy)
    }
  }

  cbind(mg_grp[, c("cell_id", "x", "y", "age", "replicate", "sig_norm")], as.data.frame(res))
}

# Split ReC by age/replicate and compute distances in parallel
mg_split <- split(mg, interaction(mg$age, mg$replicate, drop = TRUE))
dist_list <- mclapply(mg_split, compute_min_dists_for_group, all_cells = md, target_types = target_types, mc.cores = cores)
all_distances <- as.data.frame(data.table::rbindlist(dist_list))
rownames(all_distances) <- all_distances$cell_id

saveRDS(all_distances, file = file.path(out_dir, "ReC_min_distances.rds"))

# ---------- Bin signature vs. distance ----------
bin_curve <- function(dist_vec, score_vec, max_d, step, win) {
  # sliding window; returns data.frame(distance, mean, sem, n)
  mids <- numeric()
  means <- numeric()
  sems <- numeric()
  ns <- integer()
  for (d0 in seq(0, max_d - win, by = step)) {
    in_bin <- which(dist_vec >= d0 & dist_vec < (d0 + win))
    n <- length(in_bin)
    if (n > 0) {
      mids <- c(mids, d0 + win / 2)
      m <- mean(score_vec[in_bin], na.rm = TRUE)
      s <- sd(score_vec[in_bin], na.rm = TRUE)
      means <- c(means, m)
      sems <- c(sems, if (n > 1) s / sqrt(n) else 0)
      ns <- c(ns, n)
    }
  }
  data.frame(distance = mids, mean = means, sem = sems, n = ns)
}

# Build curves per (age, target_type)
curves <- list()
for (ag in sort(unique(all_distances$age))) {
  df_ag <- all_distances %>% filter(age == ag)
  for (tt in target_types) {
    dcol <- paste0("dist_to_", tt)
    keep <- !is.na(df_ag[[dcol]])
    if (!any(keep)) next
    curve <- bin_curve(
      dist_vec = df_ag[[dcol]][keep],
      score_vec = df_ag$sig_norm[keep],
      max_d = bin_max, step = bin_step, win = bin_window
    )
    if (nrow(curve) > 0) {
      curve$age <- ag
      curve$cell_type <- tt
      curves[[paste(ag, tt, sep = "_")]] <- curve
    }
  }
}
curves_df <- as.data.frame(data.table::rbindlist(curves, fill = TRUE))

# Center per (age, cell_type) so y=0 means “overall mean”
curves_df <- curves_df %>%
  group_by(age, cell_type) %>%
  mutate(mean_centered = mean - mean(mean, na.rm = TRUE)) %>%
  ungroup()

saveRDS(curves_df, file = file.path(out_dir, "signature_vs_distance_binned.rds"))

# ---------- Plots ----------
p1 <- ggplot(curves_df, aes(distance, mean_centered, color = age, fill = age)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_centered - sem, ymax = mean_centered + sem), alpha = 0.15, linewidth = 0) +
  facet_wrap(~cell_type, scales = "free_y") +
  labs(x = "Distance to nearest cell (µm)", y = "Signature (mean-centered)", title = signature_col) +
  theme_classic()

ggsave(file.path(out_dir, "signature_vs_distance_by_type.pdf"), p1, width = 9, height = 5)

# Optional: show a subset of cell types
# keep_types <- c("31 OPC-Oligo","30 Astro-Epen","29 CB Glut","28 CB GABA","33 Vascular")
# p_sub <- p1 + facet_wrap(~ cell_type, scales = "free_y", ncol = 3, labeller = as_labeller(setNames(keep_types, keep_types)))
# ggsave(file.path(out_dir, "signature_vs_distance_subset.pdf"), p_sub, width = 7, height = 4)

message("Done. Outputs written to: ", out_dir)
