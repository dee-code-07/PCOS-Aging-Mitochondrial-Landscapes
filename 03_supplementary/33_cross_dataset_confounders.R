suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(tibble)
  library(patchwork)
})

set.seed(42)

# Create output dir
out_dir <- "/home/deekshah/mini_project/analysis/21_cross_dataset_confounders/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load annotated Seurat objects
print("Loading annotated objects...")
pcos_annot <- readRDS("/home/deekshah/mini_project/scrna/output/07_annotation/pcos/pcos_annotated.rds")
aging_annot <- readRDS("/home/deekshah/mini_project/scrna/output/07_annotation/aging/aging_annotated.rds")

# 2. Load module-scored versions
print("Loading module-scored objects...")
pcos_mod <- readRDS("/home/deekshah/mini_project/analysis/03_modules/pcos_with_modules.rds")
aging_mod <- readRDS("/home/deekshah/mini_project/analysis/03_modules/aging_with_modules.rds")

gc()

# 3. Extract technical covariates and module scores
print("Extracting metadata...")
extract_meta <- function(seurat_annot, seurat_mod, dataset_name) {
  meta1 <- seurat_annot@meta.data
  meta2 <- seurat_mod@meta.data
  
  # Find doublet column
  dbl_cols <- c('scDblFinder.score', 'doublet_score', 'DF.classifications')
  dbl_col <- intersect(dbl_cols, colnames(meta1))[1]
  
  cols_to_keep <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
  if (!is.na(dbl_col)) {
    cols_to_keep <- c(cols_to_keep, dbl_col)
  }
  
  meta_tech <- meta1[, intersect(cols_to_keep, colnames(meta1)), drop = FALSE]
  if (!is.na(dbl_col) && dbl_col %in% colnames(meta_tech)) {
    colnames(meta_tech)[colnames(meta_tech) == dbl_col] <- "Doublet_Score"
    if(is.factor(meta_tech$Doublet_Score) || is.character(meta_tech$Doublet_Score)) {
        if(any(grepl("^[0-9.]+$", meta_tech$Doublet_Score))) {
            meta_tech$Doublet_Score <- as.numeric(as.character(meta_tech$Doublet_Score))
        }
    }
  }
  
  meta_tech$cell <- rownames(meta_tech)
  meta_tech$Dataset <- dataset_name
  
  # Module scores
  mod_cols <- c("Mitochondrial_Score1", "Senescence_Score1", "Oxidative_Score1", "Inflammation_Score1")
  meta_mod <- meta2[, intersect(mod_cols, colnames(meta2)), drop = FALSE]
  meta_mod$cell <- rownames(meta_mod)
  
  # merge
  merged_meta <- merge(meta_tech, meta_mod, by = "cell")
  return(merged_meta)
}

meta_pcos <- extract_meta(pcos_annot, pcos_mod, "PCOS")
meta_aging <- extract_meta(aging_annot, aging_mod, "Aging")

all_meta <- bind_rows(meta_pcos, meta_aging)
rm(pcos_annot, aging_annot, pcos_mod, aging_mod)
gc()

# 4. Compute per-dataset summary statistics
print("Computing summary statistics...")
covariates <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
if ("Doublet_Score" %in% colnames(all_meta) && is.numeric(all_meta$Doublet_Score)) {
  covariates <- c(covariates, "Doublet_Score")
}

summary_list <- lapply(covariates, function(covar) {
  if (covar %in% colnames(all_meta)) {
    all_meta %>%
      group_by(Dataset) %>%
      summarise(
        Covariate = covar,
        n_cells = n(),
        median = median(!!sym(covar), na.rm = TRUE),
        mean = mean(!!sym(covar), na.rm = TRUE),
        Q1 = quantile(!!sym(covar), 0.25, na.rm = TRUE),
        Q3 = quantile(!!sym(covar), 0.75, na.rm = TRUE),
        min = min(!!sym(covar), na.rm = TRUE),
        max = max(!!sym(covar), na.rm = TRUE)
      )
  }
})
summary_stats <- bind_rows(summary_list)

# Wilcoxon
wilcox_res <- lapply(covariates, function(covar) {
  if (covar %in% colnames(all_meta)) {
    val_pcos <- all_meta %>% filter(Dataset == "PCOS") %>% pull(!!sym(covar))
    val_aging <- all_meta %>% filter(Dataset == "Aging") %>% pull(!!sym(covar))
    w <- wilcox.test(val_pcos, val_aging)
    data.frame(Covariate = covar, Wilcox_p_value = w$p.value)
  }
}) %>% bind_rows()

summary_stats <- left_join(summary_stats, wilcox_res, by = "Covariate")

write.csv(summary_stats, file.path(out_dir, "technical_covariate_summary.csv"), row.names = FALSE)

# 5. Spearman correlations
print("Computing Spearman correlations...")
mod_cols <- c("Mitochondrial_Score1", "Senescence_Score1", "Oxidative_Score1", "Inflammation_Score1")

compute_cor <- function(df) {
  if (nrow(df) > 10000) {
    df <- df[sample(nrow(df), 10000), ]
  }
  
  cor_res <- expand.grid(Covariate = covariates, Module = mod_cols, stringsAsFactors = FALSE)
  cor_res$rho <- NA
  cor_res$p_value <- NA
  
  for (i in 1:nrow(cor_res)) {
    covar <- cor_res$Covariate[i]
    mod <- cor_res$Module[i]
    if (covar %in% colnames(df) && mod %in% colnames(df)) {
      if (is.numeric(df[[covar]]) && is.numeric(df[[mod]])) {
        test <- cor.test(df[[covar]], df[[mod]], method = "spearman", exact = FALSE)
        cor_res$rho[i] <- test$estimate
        cor_res$p_value[i] <- test$p.value
      }
    }
  }
  return(cor_res)
}

cor_pcos <- compute_cor(meta_pcos)
cor_pcos$Dataset <- "PCOS"

cor_aging <- compute_cor(meta_aging)
cor_aging$Dataset <- "Aging"

all_cor <- bind_rows(cor_pcos, cor_aging) %>%
  filter(!is.na(rho))

write.csv(all_cor, file.path(out_dir, "covariate_module_correlations.csv"), row.names = FALSE)

# 6. Generate heatmaps
print("Generating heatmaps...")
all_cor$signif <- ifelse(all_cor$p_value < 0.001, "***",
                         ifelse(all_cor$p_value < 0.01, "**",
                                ifelse(all_cor$p_value < 0.05, "*", "")))

make_heatmap <- function(cor_df, title) {
  ggplot(cor_df, aes(x = Module, y = Covariate, fill = rho)) +
    geom_tile(color = "white") +
    geom_text(aes(label = signif), color = "black", size = 5, vjust = 0.5) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, x = "", y = "")
}

p1 <- make_heatmap(all_cor %>% filter(Dataset == "PCOS"), "PCOS")
p2 <- make_heatmap(all_cor %>% filter(Dataset == "Aging"), "Aging")

combined_plot <- p1 + p2 + plot_layout(guides = "collect")

tiff(file.path(out_dir, "FigSX_covariate_correlation_heatmap.tiff"), width = 10, height = 5, units = "in", res = 600)
print(combined_plot)
dev.off()

png(file.path(out_dir, "FigSX_covariate_correlation_heatmap.png"), width = 10, height = 5, units = "in", res = 300)
print(combined_plot)
dev.off()

print("Done!")
