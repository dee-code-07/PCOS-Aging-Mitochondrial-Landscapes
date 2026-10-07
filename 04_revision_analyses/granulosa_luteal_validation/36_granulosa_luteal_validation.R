#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(ggpubr)
})

set.seed(42)

# Output directory
out_dir <- "/home/deekshah/mini_project/analysis/24_granulosa_luteal/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

###############################################################################
# Part A: Granulosa_luteal cluster composition by condition
###############################################################################
cat("Starting Part A...\n")
pcos_anno_path <- "/home/deekshah/mini_project/scrna/output/07_annotation/pcos/pcos_annotated.rds"
pcos_anno <- readRDS(pcos_anno_path)

# Add condition
pcos_anno$condition <- NA
pcos_anno$condition[grepl('case', pcos_anno$sample_id)] <- 'PCOS'
pcos_anno$condition[grepl('ctrl|control', pcos_anno$sample_id)] <- 'Control'

# Subset
gl_cells <- subset(pcos_anno, final_celltype == 'granulosa_luteal')

# Composition
total_count <- ncol(gl_cells)
cond_counts <- table(gl_cells$condition)
cond_fracs <- prop.table(cond_counts)
sample_counts <- table(gl_cells$sample_id, gl_cells$condition)

cat("Total granulosa_luteal cells:", total_count, "\n")

comp_df <- as.data.frame(table(sample = gl_cells$sample_id, condition = gl_cells$condition))
comp_df <- comp_df[comp_df$Freq > 0, ]

write.csv(comp_df, file.path(out_dir, "granulosa_luteal_composition_persample.csv"), row.names = FALSE)

# Plot
p_comp <- ggplot(comp_df, aes(x = condition, y = Freq, fill = sample)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(title = "Granulosa Luteal Cells per Sample", x = "Condition", y = "Count")

ggsave(file.path(out_dir, "granulosa_luteal_composition.png"), p_comp, width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "granulosa_luteal_composition.tiff"), p_comp, width = 6, height = 5, dpi = 600)

rm(pcos_anno, gl_cells)
gc()

###############################################################################
# Part B: Senescence module score correlation with Cdkn2a and Cdkn1a
###############################################################################
cat("Starting Part B...\n")
pcos_mod_path <- "/home/deekshah/mini_project/analysis/03_modules/pcos_with_modules.rds"
pcos_mod <- readRDS(pcos_mod_path)

# Ensure condition is added
pcos_mod$condition <- NA
pcos_mod$condition[grepl('case', pcos_mod$sample_id)] <- 'PCOS'
pcos_mod$condition[grepl('ctrl|control', pcos_mod$sample_id)] <- 'Control'

# Find senescence score column
meta_cols <- colnames(pcos_mod@meta.data)
sen_col <- meta_cols[grepl("Senescence_Score1|Senescence.*1$", meta_cols, ignore.case = TRUE)][1]

cat("Senescence column:", sen_col, "\n")

gl_mod <- subset(pcos_mod, final_celltype == 'granulosa_luteal')
gl_mod_ctrl <- subset(gl_mod, condition == 'Control')

# Function to get expression and compute correlation
get_cor_and_plot <- function(seurat_obj, gene, score_col, prefix, subset_name) {
  # Determine which assay has the gene — prefer SCT, fall back to RNA
  assay_to_use <- NULL
  if (gene %in% rownames(seurat_obj[["SCT"]])) {
    assay_to_use <- "SCT"
  } else if ("RNA" %in% Assays(seurat_obj) && gene %in% rownames(seurat_obj[["RNA"]])) {
    assay_to_use <- "RNA"
    cat("Gene", gene, "not found in SCT assay (likely excluded from top HVGs). Using RNA assay instead.\n")
  } else {
    cat("Gene", gene, "not found in SCT or RNA assay. Skipping.\n")
    return(NULL)
  }
  
  # Extract expression from the available assay
  expr <- GetAssayData(seurat_obj, assay = assay_to_use, layer = "data")[gene, ]
  scores <- seurat_obj@meta.data[[score_col]]
  
  df <- data.frame(Expression = expr, SenescenceScore = scores)
  
  res <- cor.test(df$Expression, df$SenescenceScore, method = "spearman")
  
  p <- ggscatter(df, x = "SenescenceScore", y = "Expression", 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "spearman",
                 title = paste0(gene, " vs Senescence (", subset_name, ")"),
                 xlab = "Senescence Module Score", ylab = paste0(gene, " ", assay_to_use, " Expression"))
  
  ggsave(file.path(out_dir, paste0(prefix, "_", gene, "_scatter.png")), p, width = 5, height = 5, dpi = 300)
  ggsave(file.path(out_dir, paste0(prefix, "_", gene, "_scatter.tiff")), p, width = 5, height = 5, dpi = 600)
  
  return(data.frame(Gene = gene, Subset = subset_name, Rho = res$estimate, PValue = res$p.value))
}

cor_res <- list()
for (g in c("Cdkn2a", "Cdkn1a")) {
  res_ctrl <- get_cor_and_plot(gl_mod_ctrl, g, sen_col, "Control", "Control")
  res_all <- get_cor_and_plot(gl_mod, g, sen_col, "All", "All")
  if (!is.null(res_ctrl)) cor_res[[length(cor_res) + 1]] <- res_ctrl
  if (!is.null(res_all)) cor_res[[length(cor_res) + 1]] <- res_all
}

if (length(cor_res) > 0) {
  cor_df <- do.call(rbind, cor_res)
  write.csv(cor_df, file.path(out_dir, "senescence_correlation_results.csv"), row.names = FALSE)
}

rm(gl_mod_ctrl, gl_mod)
gc()

###############################################################################
# Part C: Median module scores and proportion of cells with positive scores
###############################################################################
cat("Starting Part C...\n")

# Identify module columns
mod_cols <- meta_cols[grepl("1$", meta_cols)]
target_mods <- c("Mitochondrial", "Senescence", "Oxidative", "Inflammation")
actual_mod_cols <- sapply(target_mods, function(m) {
  mod_cols[grepl(m, mod_cols, ignore.case = TRUE)][1]
})
actual_mod_cols <- actual_mod_cols[!is.na(actual_mod_cols)]

cat("Modules found:", paste(actual_mod_cols, collapse = ", "), "\n")

calc_stats <- function(seurat_obj, cond_col) {
  meta <- seurat_obj@meta.data
  
  res_list <- lapply(actual_mod_cols, function(mc) {
    meta %>%
      group_by(final_celltype, !!sym(cond_col)) %>%
      summarise(
        Module = mc,
        MedianScore = median(!!sym(mc), na.rm = TRUE),
        MeanScore = mean(!!sym(mc), na.rm = TRUE),
        PropPositive = mean(!!sym(mc) > 0, na.rm = TRUE),
        .groups = "drop"
      )
  })
  
  do.call(rbind, res_list)
}

pcos_stats <- calc_stats(pcos_mod, "condition")
write.csv(pcos_stats, file.path(out_dir, "pcos_module_stats.csv"), row.names = FALSE)

rm(pcos_mod)
gc()

# Aging
aging_mod_path <- "/home/deekshah/mini_project/analysis/03_modules/aging_with_modules.rds"
if (file.exists(aging_mod_path)) {
  aging_mod <- readRDS(aging_mod_path)
  
  aging_mod$condition <- NA
  aging_mod$condition[grepl('_a_|aged', aging_mod$sample_id, ignore.case = TRUE)] <- 'Aged'
  aging_mod$condition[grepl('_y_|young', aging_mod$sample_id, ignore.case = TRUE)] <- 'Young'
  
  aging_stats <- calc_stats(aging_mod, "condition")
  write.csv(aging_stats, file.path(out_dir, "aging_module_stats.csv"), row.names = FALSE)
  
  rm(aging_mod)
  gc()
} else {
  cat("Aging modules file not found:", aging_mod_path, "\n")
}

cat("Done!\n")
