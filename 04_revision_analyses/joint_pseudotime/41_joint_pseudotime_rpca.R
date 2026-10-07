# ------------------------------------------------------------------------------
# Script: 29_joint_pseudotime_rpca.R
# Description: Joint pseudotime analysis using RPCA vs Harmony.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(slingshot)
  library(SingleCellExperiment)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(FNN)
})

# Increase allowed object size for Seurat parallelization
options(future.globals.maxSize = 10000 * 1024^2)
set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
is_smoke_test <- "--smoke" %in% args

# --- Configuration ---
pcos_path <- "/home/deekshah/mini_project/scrna/output/07_annotation/pcos/pcos_annotated.rds"
aging_path <- "/home/deekshah/mini_project/scrna/output/07_annotation/aging/aging_annotated.rds"
out_dir <- "/home/deekshah/mini_project/analysis/29_joint_pseudotime_rpca/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "progress_29.log")

log_msg <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
  cat(paste0("[", Sys.time(), "] ", msg, "\n"), file = log_file, append = TRUE)
}

# --- 1. Load Data ---
log_msg("Loading datasets...")
pcos_obj <- readRDS(pcos_path)
aging_obj <- readRDS(aging_path)

# --- 2. Add Metadata ---
log_msg("Formatting metadata...")
pcos_obj$dataset <- "PCOS"
pcos_obj$condition_group <- ifelse(grepl("case", pcos_obj$sample_id, ignore.case = TRUE), "PCOS_case", 
                                   ifelse(grepl("ctrl|control", pcos_obj$sample_id, ignore.case = TRUE), "PCOS_control", NA))

aging_obj$dataset <- "Aging"
aging_obj$condition_group <- ifelse(grepl("_a_|aged", aging_obj$sample_id, ignore.case = TRUE), "Aging_aged",
                                    ifelse(grepl("_y_|young", aging_obj$sample_id, ignore.case = TRUE), "Aging_young", NA))

# --- 3. Subset to Granulosa Cells ---
log_msg("Subsetting to granulosa cells...")
target_cells <- c("granulosa", "granulosa_antral", "granulosa_luteal", "granulosa_preantral")

pcos_granulosa <- subset(pcos_obj, final_celltype %in% target_cells)
aging_granulosa <- subset(aging_obj, final_celltype %in% target_cells)

rm(pcos_obj, aging_obj); gc()

# --- Smoke Test Subsetting ---
if (is_smoke_test) {
  log_msg("SMOKE TEST: Subsetting to ~3000 cells total")
  set.seed(42)
  p_cells <- sample(colnames(pcos_granulosa), min(1500, ncol(pcos_granulosa)))
  a_cells <- sample(colnames(aging_granulosa), min(1500, ncol(aging_granulosa)))
  pcos_granulosa <- subset(pcos_granulosa, cells = p_cells)
  aging_granulosa <- subset(aging_granulosa, cells = a_cells)
}

# --- 4. Merge Datasets ---
log_msg("Merging datasets...")
DefaultAssay(pcos_granulosa) <- "RNA"
DefaultAssay(aging_granulosa) <- "RNA"

for (a in setdiff(names(pcos_granulosa@assays), "RNA")) pcos_granulosa[[a]] <- NULL
for (a in setdiff(names(aging_granulosa@assays), "RNA")) aging_granulosa[[a]] <- NULL

pcos_granulosa[["RNA"]] <- JoinLayers(pcos_granulosa[["RNA"]])
aging_granulosa[["RNA"]] <- JoinLayers(aging_granulosa[["RNA"]])

merged_obj <- merge(pcos_granulosa, y = aging_granulosa)
rm(pcos_granulosa, aging_granulosa); gc()

merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])
merged_obj$sample_id <- as.character(merged_obj$sample_id)
merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], f = merged_obj$sample_id)

# --- 5. Processing & RPCA Integration ---
log_msg("Processing merged dataset (Normalize, FindVar, Scale, PCA)...")
merged_obj <- NormalizeData(merged_obj, verbose = FALSE)
merged_obj <- FindVariableFeatures(merged_obj, nfeatures = 3000, verbose = FALSE)
merged_obj <- ScaleData(merged_obj, verbose = FALSE)
merged_obj <- RunPCA(merged_obj, npcs = 30, verbose = FALSE)

log_msg("Running RPCA Integration...")
min_cells_in_sample <- min(table(merged_obj$sample_id))
k_weight <- min(30, min_cells_in_sample - 1) 
k_filter <- min(200, min_cells_in_sample - 1)

merged_obj <- IntegrateLayers(
  object = merged_obj, 
  method = RPCAIntegration, 
  orig.reduction = "pca", 
  new.reduction = "integrated.rpca",
  k.weight = k_weight,
  k.filter = k_filter,
  verbose = FALSE
)

# --- 6. Harmony Integration (for comparison) ---
log_msg("Running Harmony Integration...")
merged_obj <- RunHarmony(merged_obj, group.by.vars = "sample_id", reduction.use = "pca", reduction.save = "harmony", assay.use = "RNA", verbose = FALSE)

# --- 7. UMAP & Slingshot ---
log_msg("Joining layers back for Slingshot...")
merged_obj[["RNA"]] <- JoinLayers(merged_obj[["RNA"]])

run_trajectory <- function(obj, reduction_name, suffix) {
  log_msg(paste0("Running UMAP and Slingshot for: ", suffix))
  obj <- RunUMAP(obj, reduction = reduction_name, dims = 1:20, reduction.name = paste0("umap.", suffix), verbose = FALSE)
  
  sce <- as.SingleCellExperiment(obj, assay = "RNA")
  start_clus <- "granulosa_preantral"
  if (!(start_clus %in% unique(sce$final_celltype))) start_clus <- unique(sce$final_celltype)[1]
  
  sce <- slingshot(sce, clusterLabels = "final_celltype", reducedDim = toupper(paste0("UMAP.", suffix)), start.clus = start_clus)
  pt_matrix <- slingPseudotime(sce)
  best_lineage <- names(which.max(colSums(!is.na(pt_matrix))))
  obj[[paste0("pseudotime_", suffix)]] <- pt_matrix[, best_lineage]
  return(obj)
}

merged_obj <- run_trajectory(merged_obj, "integrated.rpca", "rpca")
merged_obj <- run_trajectory(merged_obj, "harmony", "harmony")

# --- 8. Concordance & Mixing Metrics ---
log_msg("Computing concordance and mixing metrics...")

# Spearman correlation of pseudotime
df <- merged_obj@meta.data
valid_cells <- df[!is.na(df$pseudotime_rpca) & !is.na(df$pseudotime_harmony), ]
spearman_cor <- cor(valid_cells$pseudotime_rpca, valid_cells$pseudotime_harmony, method = "spearman")

# Dataset mixing fraction (k=30 nearest neighbors)
compute_mixing <- function(reduction_name) {
  coords <- Embeddings(merged_obj, reduction_name)
  knn <- get.knn(coords, k = 30)
  datasets <- merged_obj$dataset
  mixing_fractions <- sapply(1:nrow(coords), function(i) {
    self_dataset <- datasets[i]
    neighbor_datasets <- datasets[knn$nn.index[i, ]]
    sum(neighbor_datasets != self_dataset) / 30
  })
  return(mean(mixing_fractions))
}

mix_rpca <- compute_mixing("integrated.rpca")
mix_harmony <- compute_mixing("harmony")

metrics_df <- data.frame(
  Metric = c("Spearman_Correlation", "Mean_Mixing_RPCA", "Mean_Mixing_Harmony"),
  Value = c(spearman_cor, mix_rpca, mix_harmony)
)
write.csv(metrics_df, file.path(out_dir, "integration_metrics.csv"), row.names = FALSE)
log_msg(paste0("Spearman Cor: ", round(spearman_cor, 4), " | Mixing RPCA: ", round(mix_rpca, 4), " | Mixing Harmony: ", round(mix_harmony, 4)))

# --- 9. Wilcoxon pairwise (RPCA) ---
log_msg("Computing pairwise Wilcoxon for RPCA...")
pw_res <- pairwise.wilcox.test(df$pseudotime_rpca, df$condition_group, p.adjust.method = "BH")
pw_df <- as.data.frame(as.table(pw_res$p.value))
colnames(pw_df) <- c("Group1", "Group2", "p.adj")
pw_df <- pw_df %>% filter(!is.na(p.adj))
write.csv(pw_df, file.path(out_dir, "pseudotime_pairwise_wilcoxon_rpca.csv"), row.names = FALSE)

# --- 10. Checkpoint ---
log_msg("Saving checkpoint...")
saveRDS(merged_obj, file.path(out_dir, "merged_rpca_harmony_checkpoint.rds"))

log_msg("✅ Script 29 completed.")
