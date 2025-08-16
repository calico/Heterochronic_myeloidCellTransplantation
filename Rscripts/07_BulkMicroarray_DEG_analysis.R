#### 0) Setup ------------------------------------------------------------------
suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

#### 1) Download & load GSE62420 (Affy HT MG-430 PM, GPL11180) ----------------
gset_list <- getGEO("GSE62420", GSEMatrix = TRUE, AnnotGPL = TRUE)
idx <- if (length(gset_list) > 1) grep("GPL11180", names(gset_list)) else 1
gset <- gset_list[[idx]]

# standardize feature names for convenience
fvarLabels(gset) <- make.names(fvarLabels(gset))

#### 2) Build groups from metadata --------------------------------------------
pd <- pData(gset)

# Try to locate region/tissue and age fields (handles common GEO naming styles)
find_col <- function(df, patterns) {
  hits <- lapply(patterns, function(p) grep(p, tolower(colnames(df))))
  unique(unlist(hits))
}
region_col_candidates <- find_col(pd, c("region", "brain", "tissue"))
age_col_candidates    <- find_col(pd, c("age", "month", "mo"))

region_raw <- if (length(region_col_candidates)) pd[[ region_col_candidates[1] ]] else NA
age_raw    <- if (length(age_col_candidates))    pd[[ age_col_candidates[1]    ]] else NA

# Clean to compact labels like: Ceb_04M, Cor_12M, ...
norm_region <- function(x) {
  x <- tolower(as.character(x))
  x <- case_when(
    str_detect(x, "cerebell") ~ "Ceb",
    str_detect(x, "cortex|cor") ~ "Cor",
    str_detect(x, "hipp|hip") ~ "Hip",
    str_detect(x, "striat|str") ~ "Str",
    TRUE ~ NA_character_
  )
  x
}
norm_age <- function(x) {
  x <- tolower(as.character(x))
  n <- str_match(x, "(\\d+)[ ]*m")[,2]  # pull number before 'm'/'mo'/'months'
  n <- ifelse(is.na(n), str_match(x, "(\\d+)")[,2], n)
  n <- ifelse(is.na(n), NA, sprintf("%02dM", as.integer(n)))
  n
}

region <- norm_region(region_raw)
age    <- norm_age(age_raw)

# Fallback if parsing failed: use the curated string from the original script
if (any(is.na(region)) || any(is.na(age))) {
  # map encoded groups “AAAABBBB…” -> 12 × groups across ages
  # (letters here are a curated mapping shipped with the paper’s re-analysis)
  gsms <- "AAAABBBBCCCCDDDDXXXXXXXXEEEEFFFFGGGGHHHHIIIIJJJJKKKKLLLL"
  sml  <- strsplit(gsms, "")[[1]]
  keep <- which(sml != "X")
  gset <- gset[, keep]
  sml  <- sml[keep]
  # order: Ceb_04M, Cor_04M, Hip_04M, Str_04M, Ceb_12M, Cor_12M, Hip_12M, Str_12M, Ceb_22M, Cor_22M, Hip_22M, Str_22M
  groups <- make.names(c("Ceb_04M","Cor_04M","Hip_04M","Str_04M",
                         "Ceb_12M","Cor_12M","Hip_12M","Str_12M",
                         "Ceb_22M","Cor_22M","Hip_22M","Str_22M"))
  gs <- factor(sml, labels = groups[seq_len(nlevels(factor(sml)))])
} else {
  gs <- factor(paste0(region, "_", age))
}

gset$group <- gs
levels(gs)
#> should include Ceb_04M, Cor_04M, Hip_04M, Str_04M, Ceb_12M, ... , Str_22M

#### 3) Log2 + between-array normalization ------------------------------------
ex <- exprs(gset)
qx <- quantile(ex, c(0, .25, .5, .75, .99, 1), na.rm = TRUE)
need_log2 <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)
if (need_log2) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}
ex <- normalizeBetweenArrays(ex)
exprs(gset) <- ex

# drop probes with any NA (simplifies limma)
keep <- complete.cases(ex)
gset <- gset[keep, ]
ex   <- exprs(gset)

#### 4) limma design & all pairwise contrasts ---------------------------------
design <- model.matrix(~ 0 + gs)      # one column per group
colnames(design) <- levels(gs)

fit <- lmFit(ex, design)
# all directional pairwise combinations
pairs <- combn(levels(gs), 2, simplify = FALSE)
contrast_terms <- unlist(lapply(pairs, function(p) c(paste(p[1], p[2], sep = "-"),
                                                     paste(p[2], p[1], sep = "-"))))
cont.matrix <- makeContrasts(contrasts = contrast_terms, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont.matrix))

# classify DE probes for each contrast (BH 0.05), no logFC threshold
dT <- decideTests(fit2, adjust.method = "fdr", p.value = 0.05)

#### 5) MA plot for young (04M) Cerebellum vs Cortex --------------------------
cmp <- "Ceb_04M-Cor_04M"
stopifnot(cmp %in% colnames(fit2$coefficients))  # (safe guard inside script ok to remove if you prefer)

ma_df <- data.frame(
  Amean   = fit2$Amean,
  logFC   = fit2$coefficients[, cmp],
  status  = factor(as.integer(dT[, cmp]), levels = c(-1, 0, 1),
                   labels = c("Down","NS","Up")),
  Gene    = fData(gset)$Gene.symbol
)

cols <- c(Down = "#377eb8", NS = "grey80", Up = "#e41a1c")

ggplot(ma_df, aes(x = Amean, y = logFC, color = status)) +
  geom_point(size = 0.7, alpha = 0.8) +
  scale_color_manual(values = cols) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
  labs(title = "MA plot: Cerebellum (04M) vs Cortex (04M)",
       x = "Mean expression (A)", y = "log2 fold-change (Ceb04M - Cor04M)",
       color = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right")
# (object is printed to the device; nothing is written to disk)
