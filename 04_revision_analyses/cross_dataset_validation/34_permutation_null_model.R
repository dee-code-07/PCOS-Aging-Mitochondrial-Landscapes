# IMPORTANT: This script is designed for HPC (>16GB RAM)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

set.seed(42)
options(future.globals.maxSize = Inf)
future::plan("sequential")

out_dir <- "/home/deekshah/mini_project/analysis/22_permutation_null/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load objects
cat("Loading datasets...\n")
pcos_obj <- readRDS("/home/deekshah/mini_project/scrna/output/07_annotation/pcos/pcos_annotated.rds")
aging_obj <- readRDS("/home/deekshah/mini_project/scrna/output/07_annotation/aging/aging_annotated.rds")

# 2. Add condition column
pcos_obj$condition <- ifelse(grepl('case', pcos_obj$sample_id), 'PCOS', 
                             ifelse(grepl('ctrl|control', pcos_obj$sample_id), 'Control', NA))
aging_obj$condition <- ifelse(grepl('_a_|aged', aging_obj$sample_id), 'Aged', 
                              ifelse(grepl('_y_|young', aging_obj$sample_id), 'Young', NA))

# Filter out NA conditions if any
pcos_obj <- subset(pcos_obj, cells = colnames(pcos_obj)[!is.na(pcos_obj$condition)])
aging_obj <- subset(aging_obj, cells = colnames(aging_obj)[!is.na(aging_obj$condition)])

# 3. Ensure RNA assay is default, JoinLayers if Assay5, NormalizeData if data layer is empty
prepare_obj <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  
  if (inherits(obj[["RNA"]], "Assay5")) {
    obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
  }
  
  has_data <- tryCatch({
    length(LayerData(obj, assay = "RNA", layer = "data")) > 0
  }, error = function(e) FALSE)
  
  if (!has_data) {
    obj <- NormalizeData(obj)
  }
  return(obj)
}

cat("Preparing objects...\n")
pcos_obj <- prepare_obj(pcos_obj)
aging_obj <- prepare_obj(aging_obj)

# 4. Subsample to max 5000 cells per condition per dataset
subsample_obj <- function(obj, condition_col) {
  cells_use <- c()
  for (cond in unique(obj[[condition_col, drop=TRUE]])) {
    cells_cond <- colnames(obj)[obj[[condition_col, drop=TRUE]] == cond]
    if (length(cells_cond) > 5000) {
      cells_cond <- sample(cells_cond, 5000)
    }
    cells_use <- c(cells_use, cells_cond)
  }
  obj <- subset(obj, cells = cells_use)
  
  obj <- FindVariableFeatures(obj, nfeatures = 5000)
  return(obj)
}

cat("Subsampling and finding variable features...\n")
pcos_sub <- subsample_obj(pcos_obj, "condition")
aging_sub <- subsample_obj(aging_obj, "condition")

rm(pcos_obj, aging_obj)
gc()

# 5. Real DEG analysis on subsampled data
cat("Running observed DEG analysis on subsampled data...\n")
Idents(pcos_sub) <- pcos_sub$condition
pcos_deg <- FindMarkers(pcos_sub, ident.1 = 'PCOS', ident.2 = 'Control', assay = "RNA", test.use = "wilcox", min.pct = 0.1)
pcos_sig <- rownames(pcos_deg[pcos_deg$p_val_adj < 0.05 & abs(pcos_deg$avg_log2FC) > 0.25, ])

Idents(aging_sub) <- aging_sub$condition
aging_deg <- FindMarkers(aging_sub, ident.1 = 'Aged', ident.2 = 'Young', assay = "RNA", test.use = "wilcox", min.pct = 0.1)
aging_sig <- rownames(aging_deg[aging_deg$p_val_adj < 0.05 & abs(aging_deg$avg_log2FC) > 0.25, ])

obs_overlap <- length(intersect(pcos_sig, aging_sig))
cat(sprintf("Observed overlap on subsampled data: %d\n", obs_overlap))

# 6. Run N_PERM = 100 permutations
N_PERM <- 100
null_overlaps <- numeric(N_PERM)

cat("Running permutations...\n")
for (i in 1:N_PERM) {
  if (i %% 10 == 0) cat(sprintf("Permutation %d of %d...\n", i, N_PERM))
  
  overlap_count <- tryCatch({
    pcos_shuffled <- sample(pcos_sub$condition)
    aging_shuffled <- sample(aging_sub$condition)
    
    Idents(pcos_sub) <- pcos_shuffled
    p_deg <- FindMarkers(pcos_sub, ident.1 = 'PCOS', ident.2 = 'Control', assay = "RNA", test.use = "wilcox", min.pct = 0.1, verbose = FALSE)
    
    Idents(aging_sub) <- aging_shuffled
    a_deg <- FindMarkers(aging_sub, ident.1 = 'Aged', ident.2 = 'Young', assay = "RNA", test.use = "wilcox", min.pct = 0.1, verbose = FALSE)
    
    p_sig <- rownames(p_deg[p_deg$p_val_adj < 0.05 & abs(p_deg$avg_log2FC) > 0.25, ])
    a_sig <- rownames(a_deg[a_deg$p_val_adj < 0.05 & abs(a_deg$avg_log2FC) > 0.25, ])
    
    length(intersect(p_sig, a_sig))
  }, error = function(e) {
    cat(sprintf("Error in permutation %d: %s\n", i, conditionMessage(e)))
    return(NA)
  })
  
  null_overlaps[i] <- overlap_count
  
  closeAllConnections()
  gc()
}

null_overlaps_clean <- na.omit(null_overlaps)

# 7. Compute statistics
mean_null <- mean(null_overlaps_clean)
sd_null <- sd(null_overlaps_clean)
if(is.na(sd_null) || sd_null == 0) sd_null <- 1e-6
empirical_p <- sum(null_overlaps_clean >= obs_overlap) / length(null_overlaps_clean)
z_score <- (obs_overlap - mean_null) / sd_null

cat("\n--- Final Results ---\n")
cat(sprintf("Observed overlap: %d\n", obs_overlap))
cat(sprintf("Mean null: %.2f\n", mean_null))
cat(sprintf("SD null: %.2f\n", sd_null))
cat(sprintf("Empirical p-value: %.5f\n", empirical_p))
cat(sprintf("Z-score: %.2f\n", z_score))

# 8. Generate outputs
summary_df <- data.frame(
  Metric = c("Observed Overlap", "Mean Null", "SD Null", "Empirical P-value", "Z-score"),
  Value = c(obs_overlap, mean_null, sd_null, empirical_p, z_score)
)
write.csv(summary_df, file.path(out_dir, "permutation_results.csv"), row.names = FALSE)

write.csv(data.frame(Permutation = 1:N_PERM, Overlap = null_overlaps), 
          file.path(out_dir, "permutation_null_distribution.csv"), row.names = FALSE)

p <- ggplot(data.frame(Overlap = null_overlaps_clean), aes(x = Overlap)) +
  geom_histogram(binwidth = 1, fill = "gray80", color = "black") +
  geom_vline(xintercept = obs_overlap, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = obs_overlap, y = Inf, label = sprintf("Obs=%d\nEmp P=%.4f", obs_overlap, empirical_p),
           vjust = 1.5, hjust = -0.1, color = "red") +
  theme_minimal() +
  labs(title = "Permutation Null Distribution of Shared DEGs",
       x = "Number of Shared DEGs", y = "Frequency")

ggsave(file.path(out_dir, "permutation_null_histogram.png"), p, width = 6, height = 5, dpi = 300)
ggsave(file.path(out_dir, "permutation_null_histogram.tiff"), p, width = 6, height = 5, dpi = 600)

cat("All permutation null model outputs saved to", out_dir, "\n")
