suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(future)
})

set.seed(42)
options(future.globals.maxSize = Inf)
plan("sequential")

out_dir <- "/home/deekshah/mini_project/analysis/23_rna_sct_sensitivity/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
obj_path <- "/home/deekshah/mini_project/scrna/output/07_annotation/pcos/pcos_annotated.rds"
if (!file.exists(obj_path)) {
  stop(paste("File not found:", obj_path))
}
cat("Loading Seurat object...\n")
seurat_obj <- readRDS(obj_path)

# Add condition column
cat("Adding condition column...\n")
seurat_obj$condition <- ifelse(grepl("case", seurat_obj$sample_id, ignore.case = TRUE), "PCOS", 
                               ifelse(grepl("ctrl|control", seurat_obj$sample_id, ignore.case = TRUE), "Control", NA))

# Ensure condition is not NA
if (any(is.na(seurat_obj$condition))) {
  warning("Some cells have NA condition. They will be excluded from DE testing.")
}

# Cell types to test
target_celltypes <- c("granulosa_antral", "stroma", "stroma_fibroblast", "phagocytes")
available_celltypes <- unique(seurat_obj$final_celltype)
celltypes_to_test <- intersect(target_celltypes, available_celltypes)

# If stroma or stroma_fibroblast is present, just pick one (the first one)
if ("stroma" %in% celltypes_to_test && "stroma_fibroblast" %in% celltypes_to_test) {
    celltypes_to_test <- celltypes_to_test[celltypes_to_test != "stroma_fibroblast"]
}
# Fallback to taking only first 3 if more than 3
if (length(celltypes_to_test) > 3) {
    celltypes_to_test <- celltypes_to_test[1:3]
}

if (length(celltypes_to_test) == 0) {
  stop("None of the specified cell types were found in final_celltype.")
}

cat("Testing cell types:", paste(celltypes_to_test, collapse = ", "), "\n")

# Seurat v5 preparation for RNA assay
if ("RNA" %in% Assays(seurat_obj) && inherits(seurat_obj[["RNA"]], "Assay5")) {
  seurat_obj <- JoinLayers(seurat_obj, assay = "RNA")
}

summary_list <- list()
gene_level_list <- list()

for (ct in celltypes_to_test) {
  cat("\nProcessing cell type:", ct, "\n")
  
  # Subset for cell type
  obj_sub <- subset(seurat_obj, subset = final_celltype == ct & !is.na(condition))
  
  # Set up Idents
  Idents(obj_sub) <- "condition"
  
  # --- RNA Assay ---
  cat("  Running RNA assay DE...\n")
  DefaultAssay(obj_sub) <- "RNA"
  # Check if data layer exists, if not normalize
  if (is.null(obj_sub[["RNA"]]$data) && inherits(obj_sub[["RNA"]], "Assay5") || 
      (!inherits(obj_sub[["RNA"]], "Assay5") && is.null(GetAssayData(obj_sub, assay = "RNA", slot = "data")))) {
      obj_sub <- NormalizeData(obj_sub, verbose = FALSE)
  }
  
  rna_res <- NULL
  tryCatch({
    rna_res <- FindMarkers(obj_sub, ident.1 = "PCOS", ident.2 = "Control", test.use = "wilcox", assay = "RNA")
    rna_res$gene <- rownames(rna_res)
    rna_res$assay <- "RNA"
  }, error = function(e) {
    cat("  Error in RNA FindMarkers:", e$message, "\n")
  })
  
  # --- SCT Assay ---
  cat("  Running SCT assay DE...\n")
  sct_res <- NULL
  
  if ("SCT" %in% Assays(obj_sub)) {
    DefaultAssay(obj_sub) <- "SCT"
    
    # Try with PrepSCTFindMarkers
    tryCatch({
      obj_sub_sct <- PrepSCTFindMarkers(obj_sub)
      sct_res <- FindMarkers(obj_sub_sct, ident.1 = "PCOS", ident.2 = "Control", test.use = "wilcox", assay = "SCT")
    }, error = function(e) {
      cat("  PrepSCTFindMarkers/FindMarkers failed:", e$message, "\n")
      cat("  Attempting downsampling for SCT...\n")
      
      tryCatch({
          # Downsample to max 5000 per condition
          cells_pcos <- which(Idents(obj_sub) == "PCOS")
          cells_ctrl <- which(Idents(obj_sub) == "Control")
          
          keep_pcos <- if (length(cells_pcos) > 5000) sample(cells_pcos, 5000) else cells_pcos
          keep_ctrl <- if (length(cells_ctrl) > 5000) sample(cells_ctrl, 5000) else cells_ctrl
          
          obj_sub_ds <- subset(obj_sub, cells = c(colnames(obj_sub)[keep_pcos], colnames(obj_sub)[keep_ctrl]))
          
          obj_sub_ds <- PrepSCTFindMarkers(obj_sub_ds)
          sct_res <<- FindMarkers(obj_sub_ds, ident.1 = "PCOS", ident.2 = "Control", test.use = "wilcox", assay = "SCT")
      }, error = function(e2) {
          cat("  Downsampled SCT DE also failed:", e2$message, "\n")
          cat("  RNA was the only viable approach for this cell type.\n")
      })
    })
    
    if (!is.null(sct_res)) {
      sct_res$gene <- rownames(sct_res)
      sct_res$assay <- "SCT"
    }
  } else {
    cat("  SCT assay not found in object.\n")
  }
  
  # --- Comparison ---
  if (!is.null(rna_res) && !is.null(sct_res)) {
    rna_degs <- rna_res %>% filter(p_val_adj < 0.05) %>% pull(gene)
    sct_degs <- sct_res %>% filter(p_val_adj < 0.05) %>% pull(gene)
    
    n_rna_degs <- length(rna_degs)
    n_sct_degs <- length(sct_degs)
    
    # Intersection and Jaccard
    intersect_degs <- intersect(rna_degs, sct_degs)
    union_degs <- union(rna_degs, sct_degs)
    jaccard <- if(length(union_degs) > 0) length(intersect_degs) / length(union_degs) else 0
    
    # Merge for correlation
    # We compare genes tested in both
    merged_res <- inner_join(
      rna_res %>% dplyr::select(gene, avg_log2FC_RNA = avg_log2FC, padj_RNA = p_val_adj),
      sct_res %>% dplyr::select(gene, avg_log2FC_SCT = avg_log2FC, padj_SCT = p_val_adj),
      by = "gene"
    )
    
    spearman_rho <- NA
    if (nrow(merged_res) > 2) {
        spearman_rho <- cor(merged_res$avg_log2FC_RNA, merged_res$avg_log2FC_SCT, method = "spearman", use = "complete.obs")
    }
    
    # Concordance of direction for shared DEGs
    if (length(intersect_degs) > 0) {
      shared_df <- merged_res %>% filter(gene %in% intersect_degs)
      concordant <- sum(sign(shared_df$avg_log2FC_RNA) == sign(shared_df$avg_log2FC_SCT), na.rm = TRUE)
      concordance_rate <- concordant / length(intersect_degs)
    } else {
      concordance_rate <- NA
    }
    
    # Summary
    summary_list[[ct]] <- data.frame(
      cell_type = ct,
      n_RNA_DEGs = n_rna_degs,
      n_SCT_DEGs = n_sct_degs,
      n_Shared_DEGs = length(intersect_degs),
      Jaccard_Index = jaccard,
      Spearman_Rho = spearman_rho,
      Direction_Concordance = concordance_rate
    )
    
    # Gene level
    merged_res$cell_type <- ct
    gene_level_list[[ct]] <- merged_res
    
  } else {
    cat("  Skipping comparison for", ct, "due to missing results.\n")
    summary_list[[ct]] <- data.frame(
      cell_type = ct,
      n_RNA_DEGs = ifelse(is.null(rna_res), NA, sum(rna_res$p_val_adj < 0.05, na.rm=TRUE)),
      n_SCT_DEGs = ifelse(is.null(sct_res), NA, sum(sct_res$p_val_adj < 0.05, na.rm=TRUE)),
      n_Shared_DEGs = NA,
      Jaccard_Index = NA,
      Spearman_Rho = NA,
      Direction_Concordance = NA
    )
  }
}

# --- Output generation ---
if (length(summary_list) > 0) {
  summary_df <- bind_rows(summary_list)
  write.csv(summary_df, file.path(out_dir, "RNA_vs_SCT_DEG_comparison.csv"), row.names = FALSE)
  cat("Saved RNA_vs_SCT_DEG_comparison.csv\n")
}

if (length(gene_level_list) > 0) {
  gene_df <- bind_rows(gene_level_list)
  write.csv(gene_df, file.path(out_dir, "RNA_vs_SCT_gene_level.csv"), row.names = FALSE)
  cat("Saved RNA_vs_SCT_gene_level.csv\n")
  
  # Plotting
  p_list <- list()
  for (ct in unique(gene_df$cell_type)) {
    df_ct <- gene_df %>% filter(cell_type == ct)
    rho <- summary_df %>% filter(cell_type == ct) %>% pull(Spearman_Rho)
    rho_text <- ifelse(is.na(rho), "NA", as.character(round(rho, 3)))
    
    p <- ggplot(df_ct, aes(x = avg_log2FC_RNA, y = avg_log2FC_SCT)) +
      geom_point(alpha = 0.5, size = 1) +
      geom_smooth(method = "lm", color = "red", se = FALSE, linetype = "dashed") +
      theme_minimal() +
      labs(title = paste0(ct),
           subtitle = paste("Spearman rho =", rho_text),
           x = "RNA avg_log2FC",
           y = "SCT avg_log2FC") +
      theme(plot.title = element_text(face = "bold"))
    
    p_list[[ct]] <- p
  }
  
  # Combine plots
  combined_plot <- wrap_plots(p_list, ncol = min(3, length(p_list)))
  
  # Save TIFF 600 DPI
  tiff_path <- file.path(out_dir, "logFC_correlation_scatter.tiff")
  ggsave(tiff_path, combined_plot, width = 12, height = 4, units = "in", dpi = 600, compression = "lzw", bg = "white")
  
  # Save PNG 300 DPI
  png_path <- file.path(out_dir, "logFC_correlation_scatter.png")
  ggsave(png_path, combined_plot, width = 12, height = 4, units = "in", dpi = 300, bg = "white")
  
  cat("Saved logFC correlation plots (TIFF & PNG).\n")
}

cat("Analysis complete.\n")
