suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(VISION)
  library(ggplot2)
  library(lsmeans)
})

set.seed(42)

project_dir <- getwd()
out_dir     <- file.path(project_dir, "outputs_scRNA")

#-------------------- CONFIG --------------------
proj_dir        <- "/PATH/TO/PROJECT"
in_rds  <- file.path(proj_dir, "input_data",
                     "_Harmony_neighborhood_withRegionLabels_mapmycells_annotated_spatial_allcluster_annotated_CerebellarLayer_annotated.rds")
out_dir        <- "/PATH/TO/output"
dir.create(out_dir, spatial_allursive = TRUE, showWarnings = FALSE)

spatial_all <- readRDS(in_rds)

spatial_all$age <- factor(spatial_all$age, levels = c("Young","Old"), ordered = TRUE)

# signature created in the 04_scRNAseq_age-related_DGE_analysis_CAAS_generation.R script
#   sig_cerebellum_shared <- make_signature(cb_up, cb_down, "Aging_Cerebellum_shared_WT_spatial_all")
sig_name <- sig_cerebellum_shared@name

#### 1) VISION scoring ---------------------------------------------------------
score_with_vision <- function(obj, sig) {
  obj[["RNA_forVISION"]] <- CreateAssayObject(GetAssayData(obj, assay = "RNA", slot = "counts"))
  DefaultAssay(obj) <- "RNA_forVISION"
  vobj <- Vision(obj, signatures = list(sig), assay = "RNA_forVISION",
                 projection_methods = NULL, pool = FALSE)
  vobj <- analyze(vobj)
  sc <- vobj@SigScores[, sig@name]
  obj@meta.data[[sig@name]] <- as.numeric(sc[rownames(obj@meta.data)])
  obj
}

spatial_all <- score_with_vision(spatial_all, sig_cerebellum_shared)


saveRDS(spatial_all, file.path(out_dir, "_Harmony_neighborhood_withRegionLabels_mapmycells_annotated_ReCcluster_annotated_CerebellarLayer_annotated_scScore_quantified.rds"))
