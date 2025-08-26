#### 0) Setup ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(dplyr)
  library(ggplot2)
})

project_dir <- getwd()
out_dir <- file.path(project_dir, "outputs_scRNA")
wt_path <- file.path(out_dir, "WTB6_sc_object_annotated_clean.rds")


# expected metadata: wt$age in {Young, Old}, wt$cell_class with Microglia, BAM, MacrophagePeripheral, Granulocyte, Tcell, NKcell, Ependymal, etc.
wt$age <- factor(wt$age, levels = c("Young", "Old"), ordered = TRUE)

#### 1) Split by age & create CellChat objects ---------------------------------
seurat_young <- subset(wt, subset = age == "Young")
seurat_old <- subset(wt, subset = age == "Old")

cellchat_young <- createCellChat(object = seurat_young, group.by = "cell_class")
cellchat_old <- createCellChat(object = seurat_old, group.by = "cell_class")

# mouse LR database
cellchat_db <- CellChatDB.mouse
cellchat_young@DB <- cellchat_db
cellchat_old@DB <- cellchat_db

#### 2) CellChat pipeline (helper) ---------------------------------------------
run_cellchat <- function(cc) {
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc)
  # cc <- filterCommunication(cc, min.cells = 5)   # optional
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}

cellchat_young <- run_cellchat(cellchat_young)
cellchat_old <- run_cellchat(cellchat_old)

#### 3) Figure-style panels ----------------------------------------------------
## 5A: incoming interaction counts toward Microglia (Young vs Old)
p_cnt_young <- netVisual_circle(
  cellchat_young@net$count,
  vertex.weight = as.numeric(table(cellchat_young@idents)),
  weight.scale = TRUE, label.edge = FALSE, targets.use = "Microglia",
  title.name = "Number of interactions — Young"
)

p_cnt_old <- netVisual_circle(
  cellchat_old@net$count,
  vertex.weight = as.numeric(table(cellchat_old@idents)),
  weight.scale = TRUE, label.edge = FALSE, targets.use = "Microglia",
  title.name = "Number of interactions — Old"
)

# print(p_cnt_young); print(p_cnt_old)

## 5B–C: merge and compare Young vs Old
cellchat_merged <- mergeCellChat(list(Young = cellchat_young, Old = cellchat_old), add.names = c("Young", "Old"))

# Overall counts and weights
gg_counts <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2))
gg_weights <- compareInteractions(cellchat_merged, measure = "weight", show.legend = FALSE, group = c(1, 2))
p_compare <- gg_counts + gg_weights
print(p_compare)

# Pathway ranking toward Microglia (stacked bars by source)
p_rank_microglia <- rankNet(
  cellchat_merged,
  mode = "comparison", stacked = TRUE,
  targets.use = "Microglia",
  sources.use = c("BAM", "ChoroidPlexusEpendymal", "Granulocyte", "MacrophagePeripheral", "NKcell", "Tcell")
)
print(p_rank_microglia)

# Example pathway bubble (IFN-II)
pathway_to_show <- "IFN-II" # use any name present in cellchat_merged@netP$pathways
p_bubble_ifn <- netVisual_bubble(cellchat_merged, signaling = pathway_to_show, comparison = c(1, 2))
print(p_bubble_ifn)
