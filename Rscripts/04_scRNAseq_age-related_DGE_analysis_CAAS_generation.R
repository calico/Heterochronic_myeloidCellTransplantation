#### 0) Setup ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(VISION)
})


project_dir <- getwd()
out_dir <- file.path(project_dir, "outputs_scRNA")

# EDIT: paths to the two processed/annotated objects from previous steps
wt_path <- file.path(out_dir, "WTB6_sc_object_annotated_clean.rds")
rec_path <- file.path(out_dir, "ReCB6_sc_object_annotated_clean.rds")

wt <- readRDS(wt_path)
rec <- readRDS(rec_path)

# Ensure expected metadata exist
wt$age <- factor(wt$age, levels = c("Young", "Old"), ordered = TRUE)
rec$age <- factor(rec$age, levels = c("Young", "Old"), ordered = TRUE)

#### 1) Helpers ----------------------------------------------------------------
run_age_deg_in_tissue <- function(obj, cell_class_label, tissue_label,
                                  test_use = "MAST",
                                  min_pct = 0.05, lfc_thr = 0.05,
                                  assay_use = "RNA") {
  x <- subset(obj, subset = cell_class == cell_class_label &
    tissue == tissue_label &
    age %in% c("Young", "Old"))

  DefaultAssay(x) <- assay_use
  Idents(x) <- "age"

  tab <- FindMarkers(
    object          = x,
    ident.1         = "Old",
    ident.2         = "Young",
    test.use        = test_use,
    min.pct         = min_pct,
    logfc.threshold = lfc_thr,
    assay           = assay_use,
    verbose         = TRUE
  )

  tab <- tab %>%
    mutate(
      p_val_BH = p.adjust(p_val, method = "BH"),
      diff_pct = pct.1 - pct.2, # Old - Young
      gene     = rownames(.)
    )

  return(tab)
}

sig_filter <- function(df, padj = 0.05, lfc = 0.25) {
  df %>% filter(p_val_BH <= padj, abs(avg_log2FC) >= lfc)
}

make_signature <- function(up_genes, down_genes, name) {
  vec <- c(rep(1, length(up_genes)), rep(-1, length(down_genes)))
  names(vec) <- c(up_genes, down_genes)
  # drop duplicates just in case
  vec <- vec[!(duplicated(names(vec)) | duplicated(names(vec), fromLast = TRUE))]
  createGeneSignature(name = name, sigData = vec)
}

#### 2) DEGs: WT microglia & ReC per tissue -----------------------------------
# EDIT: the label used for reconstituted cells in your metadata

# WT microglia
deg_wt_cortex <- run_age_deg_in_tissue(wt, "Microglia", tissue_label = "Cortex")
deg_wt_cereb <- run_age_deg_in_tissue(wt, "Microglia", tissue_label = "Cerebellum")

# Reconstituted cells
deg_rec_cortex <- run_age_deg_in_tissue(rec, "ReC", tissue_label = "Cortex")
deg_rec_cereb <- run_age_deg_in_tissue(rec, "ReC", tissue_label = "Cerebellum")

# “significant” (for plotting/signatures)
wt_cx_sig <- sig_filter(deg_wt_cortex)
rec_cx_sig <- sig_filter(deg_rec_cortex)
wt_cb_sig <- sig_filter(deg_wt_cereb)
rec_cb_sig <- sig_filter(deg_rec_cereb)

#### 3) Fig. 3D/E-style scatterplots ------------------------------------------
# Cortex: merge by gene (keep padj<=0.05 on both, plot all matches; highlight sig in both)
cx_merge_all <- inner_join(
  deg_wt_cortex %>% filter(p_val_BH <= 0.05) %>% select(gene, avg_log2FC),
  deg_rec_cortex %>% filter(p_val_BH <= 0.05) %>% select(gene, avg_log2FC),
  by = "gene",
  suffix = c(".wt", ".rec")
)

cx_merge_sig <- inner_join(
  wt_cx_sig %>% select(gene, avg_log2FC),
  rec_cx_sig %>% select(gene, avg_log2FC),
  by = "gene",
  suffix = c(".wt", ".rec")
)

p_scatter_cortex <- ggplot(cx_merge_all, aes(x = avg_log2FC.wt, y = avg_log2FC.rec)) +
  geom_point(size = 0.8, alpha = 0.8) +
  geom_point(data = cx_merge_sig, color = "red", size = 0.8) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  scale_x_continuous(limits = c(-2, 4), oob = squish, expand = c(0, 0)) +
  scale_y_continuous(limits = c(-2, 4), oob = squish, expand = c(0, 0)) +
  labs(
    x = "log2FC (Old vs Young) — WT microglia (Cortex)",
    y = "log2FC (Old vs Young) — Reconstituted (Cortex)"
  ) +
  theme_classic(base_size = 11)

# Cerebellum
cb_merge_all <- inner_join(
  deg_wt_cereb %>% filter(p_val_BH <= 0.05) %>% select(gene, avg_log2FC),
  deg_rec_cereb %>% filter(p_val_BH <= 0.05) %>% select(gene, avg_log2FC),
  by = "gene",
  suffix = c(".wt", ".rec")
)

cb_merge_sig <- inner_join(
  wt_cb_sig %>% select(gene, avg_log2FC),
  rec_cb_sig %>% select(gene, avg_log2FC),
  by = "gene",
  suffix = c(".wt", ".rec")
)

p_scatter_cerebellum <- ggplot(cb_merge_all, aes(x = avg_log2FC.wt, y = avg_log2FC.rec)) +
  geom_point(size = 0.8, alpha = 0.8) +
  geom_point(data = cb_merge_sig, color = "red", size = 0.8) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  scale_x_continuous(limits = c(-2, 4), oob = squish, expand = c(0, 0)) +
  scale_y_continuous(limits = c(-2, 4), oob = squish, expand = c(0, 0)) +
  labs(
    x = "log2FC (Old vs Young) — WT microglia (Cerebellum)",
    y = "log2FC (Old vs Young) — Reconstituted (Cerebellum)"
  ) +
  theme_classic(base_size = 11)

# In RStudio: print(p_scatter_cortex); print(p_scatter_cerebellum)

#### 4) VISION signatures (cortex, cerebellum, shared) -------------------------
# Direction-consistent overlaps (Old↑ or Old↓ in BOTH backgrounds)
cx_up <- intersect(
  wt_cx_sig %>% filter(avg_log2FC >= 0.25) %>% pull(gene),
  rec_cx_sig %>% filter(avg_log2FC >= 0.25) %>% pull(gene)
)
cx_down <- intersect(
  wt_cx_sig %>% filter(avg_log2FC <= -0.25) %>% pull(gene),
  rec_cx_sig %>% filter(avg_log2FC <= -0.25) %>% pull(gene)
)

cb_up <- intersect(
  wt_cb_sig %>% filter(avg_log2FC >= 0.25) %>% pull(gene),
  rec_cb_sig %>% filter(avg_log2FC >= 0.25) %>% pull(gene)
)
cb_down <- intersect(
  wt_cb_sig %>% filter(avg_log2FC <= -0.25) %>% pull(gene),
  rec_cb_sig %>% filter(avg_log2FC <= -0.25) %>% pull(gene)
)

# Shared across regions (direction-consistent in both cortex and cerebellum)
shared_up <- intersect(cx_up, cb_up)
shared_down <- intersect(cx_down, cb_down)

sig_cortex_shared <- make_signature(cx_up, cx_down, "Aging_Cortex_shared_WT_ReC")
sig_cerebellum_shared <- make_signature(cb_up, cb_down, "Aging_Cerebellum_shared_WT_ReC")
sig_shared_both <- make_signature(shared_up, shared_down, "Aging_SharedAcrossRegions_WT_ReC")


deg_wt_cortex$comparison <- "Old vs Young"
deg_wt_cortex$tissue <- "Cortex"
deg_wt_cortex$celltype_tested <- "Microglia"

deg_wt_cereb$comparison <- "Old vs Young"
deg_wt_cereb$tissue <- "Cerebellum"
deg_wt_cereb$celltype_tested <- "Microglia"

deg_rec_cortex$comparison <- "Young-to-Old vs Young-to-Young"
deg_rec_cortex$tissue <- "Cortex"
deg_rec_cortex$celltype_tested <- "ReC"

deg_rec_cereb$comparison <- "Young-to-Old vs Young-to-Young"
deg_rec_cereb$tissue <- "Cerebellum"
deg_rec_cereb$celltype_tested <- "ReC"

output_df <- rbind(deg_wt_cortex, deg_wt_cereb, deg_rec_cortex, deg_rec_cereb)

write.table(output_df, col.names = NA, sep = "\t", file = "/group/ohahnlab/OHahn/RNA-seq/Documentation/2024_04_17_HSCpilot/_Publication_material/Supplementary_tables/Fig3_C_D_DEGs_WTaging_ReCAging.txt")
