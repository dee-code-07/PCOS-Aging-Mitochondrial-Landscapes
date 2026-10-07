# ------------------------------------------------------------------------------
# Script: 27_joint_pseudotime.R
# Description: Joint pseudotime analysis of granulosa cells from PCOS and Aging datasets.
#              Addresses Major Concern 5 (cross-dataset pseudotime normalization).
# Note: Designed for HPC execution (High RAM usage expected, >16GB).
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(slingshot)
  library(SingleCellExperiment)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

set.seed(42)

# --- Configuration ---
pcos_path <- "/home/deekshah/mini_project/scrna/output/07_annotation/pcos/pcos_annotated.rds"
aging_path <- "/home/deekshah/mini_project/scrna/output/07_annotation/aging/aging_annotated.rds"
out_dir <- "/home/deekshah/mini_project/analysis/27_joint_pseudotime/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load Data ---
cat("Loading datasets...\n")
pcos_obj <- readRDS(pcos_path)
aging_obj <- readRDS(aging_path)

# --- 2. Add Metadata ---
cat("Formatting metadata...\n")
pcos_obj$dataset <- "PCOS"
pcos_obj$condition_group <- ifelse(grepl("case", pcos_obj$sample_id, ignore.case = TRUE), "PCOS_case", 
                                   ifelse(grepl("ctrl|control", pcos_obj$sample_id, ignore.case = TRUE), "PCOS_control", NA))

aging_obj$dataset <- "Aging"
aging_obj$condition_group <- ifelse(grepl("_a_|aged", aging_obj$sample_id, ignore.case = TRUE), "Aging_aged",
                                    ifelse(grepl("_y_|young", aging_obj$sample_id, ignore.case = TRUE), "Aging_young", NA))

# --- 3. Subset to Granulosa Cells ---
cat("Subsetting to granulosa cells...\n")
target_cells <- c("granulosa", "granulosa_antral", "granulosa_luteal", "granulosa_preantral")

pcos_granulosa <- subset(pcos_obj, final_celltype %in% target_cells)
aging_granulosa <- subset(aging_obj, final_celltype %in% target_cells)

rm(pcos_obj, aging_obj)
gc()

# --- 4. Merge Datasets ---
cat("Merging datasets...\n")
# Remove ALL non-RNA assays before merge — the two datasets have different
# SCT models/features which causes rbind column mismatch during merge.
# We re-normalize from RNA anyway.
DefaultAssay(pcos_granulosa) <- "RNA"
DefaultAssay(aging_granulosa) <- "RNA"

# Use names(obj@assays) instead of Assays() — Seurat v5 Assays() returns
# an S4 object that setdiff() cannot coerce to a vector.
pcos_assays <- names(pcos_granulosa@assays)
aging_assays <- names(aging_granulosa@assays)
cat(sprintf("PCOS assays: %s\n", paste(pcos_assays, collapse = ", ")))
cat(sprintf("Aging assays: %s\n", paste(aging_assays, collapse = ", ")))

for (a in setdiff(pcos_assays, "RNA")) {
  cat(sprintf("  Removing assay '%s' from PCOS object\n", a))
  pcos_granulosa[[a]] <- NULL
}
for (a in setdiff(aging_assays, "RNA")) {
  cat(sprintf("  Removing assay '%s' from Aging object\n", a))
  aging_granulosa[[a]] <- NULL
}

# JoinLayers before merge — Seurat v5 stores per-sample layers which cause
# 'match' requires vector arguments error during merge if not joined first.
cat("Joining layers...\n")
pcos_granulosa[["RNA"]] <- JoinLayers(pcos_granulosa[["RNA"]])
aging_granulosa[["RNA"]] <- JoinLayers(aging_granulosa[["RNA"]])

cat("Merging...\n")
merged_obj <- merge(pcos_granulosa, y = aging_granulosa)

rm(pcos_granulosa, aging_granulosa)
gc()

# --- 5. Process Merged Object ---
cat("Processing merged dataset...\n")
DefaultAssay(merged_obj) <- "RNA"
merged_obj <- JoinLayers(merged_obj)
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj, nfeatures = 3000)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj, npcs = 30)

cat("Running Harmony batch correction...\n")
merged_obj <- RunHarmony(merged_obj, group.by.vars = "sample_id", reduction.use = "pca", assay.use = "RNA")

cat("Running UMAP...\n")
merged_obj <- RunUMAP(merged_obj, reduction = "harmony", dims = 1:20)

# --- 6. Slingshot Trajectory Analysis ---
cat("Running Slingshot...\n")
sce <- as.SingleCellExperiment(merged_obj, assay = "RNA")

# Set start cluster. Fallback to just "granulosa" if "granulosa_preantral" isn't present
start_clus <- "granulosa_preantral"
if (!(start_clus %in% unique(sce$final_celltype))) {
  cat("Warning: granulosa_preantral not found. Falling back to an available cluster as start point.\n")
  avail <- unique(sce$final_celltype)
  start_clus <- avail[1] 
}

sce <- slingshot(sce, clusterLabels = "final_celltype", reducedDim = "UMAP", start.clus = start_clus)

# Extract pseudotime. If multiple lineages, take lineage 1 or the one with most cells.
pt_matrix <- slingPseudotime(sce)
# To find the lineage covering the most cells, count non-NAs
lineage_counts <- colSums(!is.na(pt_matrix))
best_lineage <- names(which.max(lineage_counts))
merged_obj$pseudotime <- pt_matrix[, best_lineage]

# --- 7. Statistical Comparisons ---
cat("Computing statistics...\n")
meta_df <- merged_obj@meta.data
meta_df$cell_id <- rownames(meta_df)
meta_df <- meta_df %>%
  filter(!is.na(pseudotime), !is.na(condition_group)) %>%
  select(cell_id, condition_group, pseudotime)

# Median pseudotime
summary_stats <- meta_df %>%
  group_by(condition_group) %>%
  summarize(
    median_pseudotime = median(pseudotime, na.rm = TRUE),
    n_cells = n()
  )

write.csv(summary_stats, file.path(out_dir, "pseudotime_summary.csv"), row.names = FALSE)

# Pairwise Wilcoxon tests
cat("Running pairwise Wilcoxon tests...\n")
pw_res <- pairwise.wilcox.test(meta_df$pseudotime, meta_df$condition_group, p.adjust.method = "BH")
pw_df <- as.data.frame(as.table(pw_res$p.value))
colnames(pw_df) <- c("Group1", "Group2", "p.adj")
pw_df <- pw_df %>% filter(!is.na(p.adj))
write.csv(pw_df, file.path(out_dir, "pseudotime_pairwise_wilcoxon.csv"), row.names = FALSE)

# --- 8. Generating Figures ---
cat("Generating figures...\n")

# UMAPs
p_cond <- DimPlot(merged_obj, group.by = "condition_group", reduction = "umap") + ggtitle("Condition Group")
p_data <- DimPlot(merged_obj, group.by = "dataset", reduction = "umap") + ggtitle("Dataset")
p_time <- FeaturePlot(merged_obj, features = "pseudotime", reduction = "umap") + 
          scale_color_viridis_c() + ggtitle("Pseudotime")

# Density
p_dens <- ggplot(meta_df, aes(x = pseudotime, fill = condition_group)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  ggtitle("Pseudotime Distribution") +
  theme(legend.position = "bottom")

# Boxplot
p_box <- ggplot(meta_df, aes(x = condition_group, y = pseudotime, fill = condition_group)) +
  geom_boxplot(alpha = 0.8) +
  theme_minimal() +
  ggtitle("Pseudotime Comparisons") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine
combined_plot <- (p_cond | p_data) / (p_time | p_dens) / (p_box | plot_spacer()) + plot_layout(heights = c(1, 1, 1))

# Save
ggsave(file.path(out_dir, "joint_pseudotime_analysis.png"), combined_plot, width = 16, height = 18, dpi = 300)
ggsave(file.path(out_dir, "joint_pseudotime_analysis.tiff"), combined_plot, width = 16, height = 18, dpi = 600, device = "tiff", compression = "lzw")

cat("Analysis complete.\n")
