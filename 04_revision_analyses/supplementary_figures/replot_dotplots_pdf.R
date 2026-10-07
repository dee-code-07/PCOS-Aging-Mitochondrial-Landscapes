#!/usr/bin/env Rscript
################################################################################
# Script: replot_dotplots_pdf.R
# Purpose: Replot dotplot_aging_top30 and dotplot_pcos_top30 as PDFs
################################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(viridis)
})

# ============================================================================
# CONFIGURATION & PATHS
# ============================================================================

PROJECT_ROOT <- "/mnt/e/Documents/mini_project"

INPUT_PRIORITY     <- file.path(PROJECT_ROOT, "analysis/09_gene_prioritization/Top_Priority_Genes_Tier1.csv")
INPUT_SEURAT_PCOS  <- file.path(PROJECT_ROOT, "scrna/output/07_annotation/pcos/pcos_annotated.rds")
INPUT_SEURAT_AGING <- file.path(PROJECT_ROOT, "scrna/output/07_annotation/aging/aging_annotated.rds")

OUTPUT_DIR <- file.path(PROJECT_ROOT, "analysis/pdf_figures")
if(!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# ============================================================================
# LOAD DATA
# ============================================================================

cat("Loading prioritized genes...\n")
priority_genes <- read_csv(INPUT_PRIORITY, show_col_types = FALSE)
top_genes <- priority_genes$gene[1:50]
plot_genes <- top_genes[1:30]

# PCOS dot plot
cat("Loading PCOS Seurat object (this may take a while)...\n")
pcos_seurat <- readRDS(INPUT_SEURAT_PCOS)
ct_col_pcos <- if("final_celltype" %in% colnames(pcos_seurat@meta.data)) "final_celltype" else "seurat_clusters"

cat("Generating PCOS dot plot...\n")
valid_genes_pcos <- intersect(plot_genes, rownames(pcos_seurat))
pcos_seurat_subset <- subset(pcos_seurat, features = valid_genes_pcos)

p1 <- DotPlot(
  pcos_seurat_subset,
  features = valid_genes_pcos,
  group.by = ct_col_pcos,
  dot.scale = 8
) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  scale_color_viridis(option = "plasma") +
  labs(title = "PCOS: Top 30 Genes by Cell Type")

ggsave(
  file.path(OUTPUT_DIR, "dotplot_pcos_top30.pdf"),
  p1,
  width = 10,
  height = 12,
  device = "pdf"
)

# Free memory
rm(pcos_seurat, pcos_seurat_subset, p1)
gc()

# Aging dot plot
cat("Loading Aging Seurat object (this may take a while)...\n")
aging_seurat <- readRDS(INPUT_SEURAT_AGING)
ct_col_aging <- if("final_celltype" %in% colnames(aging_seurat@meta.data)) "final_celltype" else "seurat_clusters"

cat("Generating Aging dot plot...\n")
valid_genes_aging <- intersect(plot_genes, rownames(aging_seurat))
aging_seurat_subset <- subset(aging_seurat, features = valid_genes_aging)

p2 <- DotPlot(
  aging_seurat_subset,
  features = valid_genes_aging,
  group.by = ct_col_aging,
  dot.scale = 8
) +
  coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  scale_color_viridis(option = "plasma") +
  labs(title = "Aging: Top 30 Genes by Cell Type")

ggsave(
  file.path(OUTPUT_DIR, "dotplot_aging_top30.pdf"),
  p2,
  width = 10,
  height = 12,
  device = "pdf"
)

cat("Done! Files saved in:", OUTPUT_DIR, "\n")
