#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(patchwork)
})

set.seed(42)

# Define paths
spatial_dir <- "/home/deekshah/mini_project/spatial/output/06_label_transfer_FIXED/aging"
module_dir <- "/home/deekshah/mini_project/scrna/resources/modules"
output_dir <- "/home/deekshah/mini_project/analysis/26_spatial_sensitivity"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load all 4 aging spatial Seurat objects
sections <- c("ya_1", "ya_2", "ya_3", "ya_4")
seurat_list <- list()

for (sec in sections) {
  obj_path <- file.path(spatial_dir, sec, paste0(sec, "_with_celltypes.rds"))
  cat("Loading", obj_path, "\n")
  obj <- readRDS(obj_path)
  obj$section <- sec
  obj$group <- ifelse(sec %in% c("ya_1", "ya_2"), "Young", "Aged")
  DefaultAssay(obj) <- "RNA"
  seurat_list[[sec]] <- obj
}

# Merge objects for easier processing
merged_obj <- merge(seurat_list[[1]], y = c(seurat_list[[2]], seurat_list[[3]], seurat_list[[4]]), add.cell.ids = sections)

# 2. Load and prepare modules
module_files <- c(
  "senescence_SASP_module.csv",
  "mitochondrial_module.csv",
  "oxidative_stress_module.csv",
  "inflammation_module.csv"
)

modules <- list()
for (mf in module_files) {
  mod_path <- file.path(module_dir, mf)
  if (file.exists(mod_path)) {
    # Read the CSV. Assuming it's a single column of genes without header, or with header.
    # We use read.csv and extract the first column.
    genes <- read.csv(mod_path, header = FALSE, stringsAsFactors = FALSE)[,1]
    
    # In case there's a header like "gene", remove it if present
    if (tolower(genes[1]) == "gene" || tolower(genes[1]) == "genes") {
        genes <- genes[-1]
    }
  } else {
    warning(paste("Module file not found:", mod_path))
    next
  }
  
  # Convert human gene symbols to mouse format
  # Special handling for mitochondrial genes: MT-CO1 -> mt-Co1, MT-ND1 -> mt-Nd1
  mouse_genes <- sapply(genes, function(g) {
    if (grepl("^MT-", g, ignore.case = TRUE)) {
      # Mitochondrial gene: lowercase "mt-" prefix, then title-case the rest
      suffix <- sub("^MT-", "", g, ignore.case = TRUE)
      paste0("mt-", toupper(substr(suffix, 1, 1)), tolower(substr(suffix, 2, nchar(suffix))))
    } else {
      # Standard gene: first letter uppercase, rest lowercase
      paste0(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g))))
    }
  }, USE.NAMES = FALSE)
  
  mod_name <- gsub("_module\\.csv$", "", mf)
  modules[[mod_name]] <- mouse_genes
}

# Score with UCell — must JoinLayers first for Seurat v5 compatibility
merged_obj <- JoinLayers(merged_obj)
merged_obj <- AddModuleScore_UCell(merged_obj, features = modules, name = "_UCell", slot = "counts")

# 3. Perform sensitivity analysis
meta <- merged_obj@meta.data

# Define comparisons
comparisons <- list(
  "ya_3_vs_Young" = list(aged = "ya_3", young = c("ya_1", "ya_2")),
  "ya_4_vs_Young" = list(aged = "ya_4", young = c("ya_1", "ya_2")),
  "Combined_vs_Young" = list(aged = c("ya_3", "ya_4"), young = c("ya_1", "ya_2"))
)

results_list <- list()
mod_names <- paste0(names(modules), "_UCell")

for (comp_name in names(comparisons)) {
  comp <- comparisons[[comp_name]]
  
  for (mod in mod_names) {
    if (!mod %in% colnames(meta)) next
    
    aged_scores <- meta[meta$section %in% comp$aged, mod]
    young_scores <- meta[meta$section %in% comp$young, mod]
    
    # Wilcoxon test
    wt <- wilcox.test(aged_scores, young_scores)
    
    # Cohen's d (pooled SD)
    n1 <- length(aged_scores)
    n2 <- length(young_scores)
    
    if (n1 > 1 && n2 > 1) {
        var1 <- var(aged_scores)
        var2 <- var(young_scores)
        mean1 <- mean(aged_scores)
        mean2 <- mean(young_scores)
        
        pooled_sd <- sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
        cohens_d <- (mean1 - mean2) / pooled_sd
    } else {
        mean1 <- NA
        mean2 <- NA
        cohens_d <- NA
    }
    
    results_list[[length(results_list) + 1]] <- data.frame(
      Comparison = comp_name,
      Module = gsub("_UCell$", "", mod),
      P_value = wt$p.value,
      Cohens_d = cohens_d,
      Mean_Aged = mean1,
      Mean_Young = mean2,
      stringsAsFactors = FALSE
    )
  }
}

results_df <- do.call(rbind, results_list)

# 4. Generate outputs
# Results CSV table
write.csv(results_df, file.path(output_dir, "sensitivity_analysis_results.csv"), row.names = FALSE)

# Grouped bar plot (Cohen's d)
p_bar <- ggplot(results_df, aes(x = Comparison, y = Cohens_d, fill = Comparison)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  facet_wrap(~ Module, scales = "free_y") +
  theme_minimal() +
  labs(title = "Effect Size (Cohen's d) by Module and Comparison",
       y = "Cohen's d (Aged vs Young)", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "cohens_d_barplot.png"), p_bar, width = 10, height = 8, dpi = 300)
ggsave(file.path(output_dir, "cohens_d_barplot.tiff"), p_bar, width = 10, height = 8, dpi = 600, compression = "lzw")

# Violin plot for Senescence
if ("senescence_SASP_UCell" %in% colnames(meta)) {
  p_violin <- ggplot(meta, aes(x = section, y = senescence_SASP_UCell, fill = group)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_classic() +
    labs(title = "Senescence SASP UCell Scores by Section",
         x = "Section", y = "UCell Score") +
    scale_fill_manual(values = c("Young" = "#00BFC4", "Aged" = "#F8766D"))

  ggsave(file.path(output_dir, "senescence_violin_plot.png"), p_violin, width = 8, height = 6, dpi = 300)
  ggsave(file.path(output_dir, "senescence_violin_plot.tiff"), p_violin, width = 8, height = 6, dpi = 600, compression = "lzw")
}

cat("Analysis complete. Results saved to", output_dir, "\n")
