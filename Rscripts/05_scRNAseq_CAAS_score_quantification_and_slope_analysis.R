suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(VISION)
  library(ggplot2)
  library(lsmeans)
})

set.seed(42)

project_dir <- getwd()
out_dir <- file.path(project_dir, "outputs_scRNA")

# same two objects as before
wt_path <- file.path(out_dir, "WTB6_sc_object_annotated_clean.rds")
rec_path <- file.path(out_dir, "ReCB6_sc_object_annotated_clean.rds")

wt <- readRDS(wt_path)
rec <- readRDS(rec_path)

wt$age <- factor(wt$age, levels = c("Young", "Old"), ordered = TRUE)
rec$age <- factor(rec$age, levels = c("Young", "Old"), ordered = TRUE)

# signature created in the previous step
#   sig_cerebellum_shared <- make_signature(cb_up, cb_down, "Aging_Cerebellum_shared_WT_ReC")
sig_name <- sig_cerebellum_shared@name

#### 1) VISION scoring ---------------------------------------------------------
score_with_vision <- function(obj, sig) {
  obj[["RNA_forVISION"]] <- CreateAssayObject(GetAssayData(obj, assay = "RNA", slot = "counts"))
  DefaultAssay(obj) <- "RNA_forVISION"
  vobj <- Vision(obj,
    signatures = list(sig), assay = "RNA_forVISION",
    projection_methods = NULL, pool = FALSE
  )
  vobj <- analyze(vobj)
  sc <- vobj@SigScores[, sig@name]
  obj@meta.data[[sig@name]] <- as.numeric(sc[rownames(obj@meta.data)])
  obj
}

wt <- score_with_vision(wt, sig_cerebellum_shared)
rec <- score_with_vision(rec, sig_cerebellum_shared)

#### 2) Slopes by group (Cortex/Cerebellum × Microglia/ReC) -------------------

## WT microglia
score_table <- wt@meta.data
score_table$score <- score_table[, sig_name]
score_table_MG <- score_table %>% filter(cell_class == "Microglia")
score_table_MG$age <- as.numeric(plyr::mapvalues(score_table_MG$age, from = c("Young", "Old"), to = c(3, 21)))

# violins + LM line (Fig. 3H style)
(myplot <- ggplot(score_table_MG, aes(x = age, y = score)) +
  geom_point(color = "#bbbdbf", position = "jitter") +
  geom_violin(notch = FALSE, outlier.colour = NA, aes(fill = as.factor(age))) +
  geom_smooth(aes(x = as.numeric(age), y = score), color = "red", se = FALSE, method = "lm") +
  facet_grid(~tissue)
)

# linear model with interaction for slopes (WT microglia)
score_table_MG$tissue <- as.character(score_table_MG$tissue)
m.interaction <- lm(score ~ age * tissue, data = score_table_MG)
m.lst <- lstrends(m.interaction, "tissue", var = "age")
slope_dataframe <- as.data.frame(m.lst)
slope_pair_wiseTestResults <- as.data.frame(pairs(m.lst))

slope_dataframe$tissue <- factor(
  slope_dataframe$tissue,
  levels = as.character(slope_dataframe[order(slope_dataframe$age.trend, decreasing = TRUE), "tissue"]),
  ordered = TRUE
)

(p_slope_bar_WT <- ggplot(slope_dataframe, aes(x = tissue, y = age.trend, fill = tissue)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.9)) +
  labs(x = NULL, y = "CAAS slope") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
)

## Reconstituted cells (label auto-detect)
rec_label <- if ("Reconstituted" %in% rec$cell_class) "Reconstituted" else "ReC"

score_table <- rec@meta.data
score_table$score <- score_table[, sig_name]
score_table_ReC <- score_table %>% filter(cell_class == rec_label)
score_table_ReC$age <- as.numeric(plyr::mapvalues(score_table_ReC$age, from = c("Young", "Old"), to = c(4, 20)))

# violins + LM line (Fig. 3H style)
(myplot <- ggplot(score_table_ReC, aes(x = age, y = score)) +
  geom_point(color = "#bbbdbf", position = "jitter") +
  geom_violin(notch = FALSE, outlier.colour = NA, aes(fill = as.factor(age))) +
  geom_smooth(aes(x = as.numeric(age), y = score), color = "red", se = FALSE, method = "lm") +
  facet_grid(~tissue)
)

# linear model with interaction for slopes (ReC)
score_table_ReC$tissue <- as.character(score_table_ReC$tissue)
m.interaction <- lm(score ~ age * tissue, data = score_table_ReC)
m.lst <- lstrends(m.interaction, "tissue", var = "age")
slope_dataframe <- as.data.frame(m.lst)
slope_pair_wiseTestResults <- as.data.frame(pairs(m.lst))

slope_dataframe$tissue <- factor(
  slope_dataframe$tissue,
  levels = as.character(slope_dataframe[order(slope_dataframe$age.trend, decreasing = TRUE), "tissue"]),
  ordered = TRUE
)

(p_slope_bar_ReC <- ggplot(slope_dataframe, aes(x = tissue, y = age.trend, fill = tissue)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(0.9)) +
  labs(x = NULL, y = "CAAS slope") +
  scale_y_continuous(expand = c(0, 0)) +
  theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
)
