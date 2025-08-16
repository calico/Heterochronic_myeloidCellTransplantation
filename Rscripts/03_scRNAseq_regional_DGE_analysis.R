#### 0) Setup ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(scales)
})

project_dir <- getwd()
out_dir     <- file.path(project_dir, "outputs_scRNA")
obj_path    <- file.path(out_dir, "sc_object_annotated_clean.rds")  # from previous step

sc_object <- readRDS(obj_path)

#### 1) Subset to WT/Young ReC ------------------------------------------
# Assumes metadata columns: age (Young/Old), tissue (Cortex/Cerebellum), cell_class (ReC)
sc_object$age <- factor(sc_object$age, levels = c("Young","Old"), ordered = TRUE)

sc_young          <- subset(sc_object, subset = age == "Young")
sc_young_mg       <- subset(sc_young, subset = cell_class == "ReC")

DefaultAssay(sc_young_mg) <- "RNA"
Idents(sc_young_mg)       <- "tissue"   # contrasts will be "Cerebellum" vs "Cortex"

#### 2) DEG (MAST) -------------------------------------------------------------
test_to_use     <- "MAST"
min_pct         <- 0.05
logfc_threshold <- 0.05
assay_to_use    <- "RNA"

deg_region <- FindMarkers(
  object          = sc_young_mg,
  ident.1         = "Cerebellum",
  ident.2         = "Cortex",
  test.use        = test_to_use,
  min.pct         = min_pct,
  logfc.threshold = logfc_threshold,
  assay           = assay_to_use,
  verbose         = TRUE
)

# tidy columns
deg_region <- deg_region %>%
  mutate(
    p_val_BH = p.adjust(p_val, method = "BH"),
    diff_pct = pct.1 - pct.2,
    gene     = rownames(.)
  )

#### 3) Simple volcano plot ----------------------------------------------------
# thresholds used in the paper-style panels
lfc_cut  <- 0.25
padj_cut <- 0.05

deg_region <- deg_region %>%
  mutate(
    status = case_when(
      p_val_BH <= padj_cut & avg_log2FC >=  lfc_cut  ~ "Up",
      p_val_BH <= padj_cut & avg_log2FC <= -lfc_cut  ~ "Down",
      TRUE                                           ~ "n.s."
    )
  )

p_volcano <- ggplot(deg_region, aes(x = avg_log2FC, y = -log10(p_val_BH))) +
  geom_point(aes(color = status), size = 1, alpha = 0.9) +
  scale_color_manual(values = c(Up = "#ea5430", Down = "#6181d1", `n.s.` = "grey80")) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 3, linewidth = 0.4) +
  geom_hline(yintercept = -log10(padj_cut), linetype = 3, linewidth = 0.4) +
  scale_x_continuous(limits = c(-3, 3), oob = squish, expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, .02))) +
  labs(x = "log2FC (Cerebellum / Cortex)", y = "-log10(adj. p-value)", color = NULL) +
  theme_classic(base_size = 11)

# Optional labels for a few key genes (edit as desired)
label_genes <- c("Tmem119","P2ry12","H2-Ab1","B2m","Il1b","Ifit3","Lpl")
lab_df <- subset(deg_region, gene %in% label_genes)
if (nrow(lab_df)) {
  p_volcano <- p_volcano +
    geom_point(data = lab_df, shape = 21, stroke = 0.6, fill = NA) +
    ggrepel::geom_text_repel(data = lab_df, aes(label = gene), size = 3, max.overlaps = 50)
}

print(p_volcano)  
