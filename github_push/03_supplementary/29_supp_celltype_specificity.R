#!/usr/bin/env Rscript
# script: 14_celltype_specificity.r
# project: mitochondrial dysfunction in pcos and ovarian aging
# purpose: determine which cell types express prioritized genes
# methods: dot plots, cell-type enrichment, specificity scores
# date: 2026-02-04
# reference: skinnider et al. (2021) genome biol pmid:33436056

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)    # Added for pivot_longer/pivot_wider
  library(readr)
  library(ggplot2)
  library(pheatmap)
  library(viridis)
})

# configuration & paths

PROJECT_ROOT <- "E:/Documents/mini_project"

# updated input paths based on your file tree
INPUT_PRIORITY     <- file.path(PROJECT_ROOT, "analysis/09_gene_prioritization/Top_Priority_Genes_Tier1.csv")
INPUT_SEURAT_PCOS  <- file.path(PROJECT_ROOT, "scrna/output/07_annotation/pcos/pcos_annotated.rds")
INPUT_SEURAT_AGING <- file.path(PROJECT_ROOT, "scrna/output/07_annotation/aging/aging_annotated.rds")

OUTPUT_DIR <- file.path(PROJECT_ROOT, "analysis/14_celltype_specificity")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# load data

cat("\n==============================================================================\n")
cat("Cell-Type Specificity Analysis\n")
cat("==============================================================================\n\n")

# check inputs
if(!file.exists(INPUT_PRIORITY)) stop("Priority genes file not found: ", INPUT_PRIORITY)
if(!file.exists(INPUT_SEURAT_PCOS)) stop("PCOS Seurat object not found: ", INPUT_SEURAT_PCOS)
if(!file.exists(INPUT_SEURAT_AGING)) stop("Aging Seurat object not found: ", INPUT_SEURAT_AGING)

priority_genes <- read_csv(INPUT_PRIORITY, show_col_types = FALSE)
top_genes <- priority_genes$gene[1:50]  # Top 50 genes

cat("Loading Seurat objects...\n")
pcos_seurat <- readRDS(INPUT_SEURAT_PCOS)
aging_seurat <- readRDS(INPUT_SEURAT_AGING)

cat("Loaded prioritized genes:", length(top_genes), "\n")
cat("PCOS cells:", ncol(pcos_seurat), "\n")
cat("Aging cells:", ncol(aging_seurat), "\n\n")

# determine correct celltype column (robust check)
ct_col_pcos <- if("final_celltype" %in% colnames(pcos_seurat@meta.data)) "final_celltype" else "seurat_clusters"
ct_col_aging <- if("final_celltype" %in% colnames(aging_seurat@meta.data)) "final_celltype" else "seurat_clusters"

cat("Using grouping columns: PCOS =", ct_col_pcos, "| Aging =", ct_col_aging, "\n\n")

# function: calculate cell-type specificity score
# tau index: 0 = ubiquitous, 1 = highly specific

calculate_tau <- function(expr_vector) {
# tau specificity index (yanai et al. 2005)
  expr_max <- max(expr_vector)
  if (expr_max == 0) return(0)
  
  tau <- sum(1 - (expr_vector / expr_max)) / (length(expr_vector) - 1)
  return(tau)
}

# calculate cell-type expression profiles

cat("Calculating cell-type expression profiles...\n")

# helper to safely join layers if using seurat v5
get_avg_expr <- function(seu, genes, group_col) {
  if (grepl("^5", packageVersion("Seurat"))) seu <- JoinLayers(seu)
  AverageExpression(seu, features = genes, group.by = group_col, slot = "data")$RNA
}

# pcos dataset
pcos_expr <- get_avg_expr(pcos_seurat, top_genes, ct_col_pcos)

pcos_expr_long <- pcos_expr %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "celltype",
    values_to = "expression"
  ) %>%
  mutate(dataset = "PCOS")

# aging dataset
aging_expr <- get_avg_expr(aging_seurat, top_genes, ct_col_aging)

aging_expr_long <- aging_expr %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "celltype",
    values_to = "expression"
  ) %>%
  mutate(dataset = "Aging")

# combine
combined_expr <- bind_rows(pcos_expr_long, aging_expr_long)

# calculate tau specificity scores

cat("Calculating Tau specificity scores...\n")

tau_scores <- combined_expr %>%
  group_by(gene, dataset) %>%
  summarise(
    tau = calculate_tau(expression),
    max_celltype = celltype[which.max(expression)],
    max_expression = max(expression),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = dataset,
    values_from = c(tau, max_celltype, max_expression),
    names_sep = "_"
  )

write_csv(
  tau_scores,
  file.path(OUTPUT_DIR, "tau_specificity_scores.csv")
)

# visualization 1: dot plot (pcos)

cat("\nGenerating dot plots...\n")

# select top 30 genes for clarity
plot_genes <- top_genes[1:30]

# pcos dot plot
# subset() can sometimes fail if genes aren't found, so we filter valid genes first
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
  file.path(OUTPUT_DIR, "dotplot_pcos_top30.png"),
  p1,
  width = 10,
  height = 12,
  dpi = 300
)

# aging dot plot
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
  file.path(OUTPUT_DIR, "dotplot_aging_top30.png"),
  p2,
  width = 10,
  height = 12,
  dpi = 300
)

# visualization 2: heatmap of expression

cat("Generating expression heatmaps...\n")

# helper to check if we have enough genes for heatmap
safe_heatmap <- function(expr_mat, genes, title, fname) {
  valid <- intersect(genes, rownames(expr_mat))
  if(length(valid) < 3) {
    cat("Skipping heatmap for", title, "- not enough valid genes.\n")
    return(NULL)
  }
  
  pheatmap(
    expr_mat[valid, , drop=FALSE],
    scale = "row",
    color = viridis(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 8,
    fontsize_col = 10,
    main = title,
    filename = fname,
    width = 8,
    height = 10
  )
}

# pcos heatmap
safe_heatmap(pcos_expr, plot_genes, 
             "PCOS: Gene Expression by Cell Type (Z-score)", 
             file.path(OUTPUT_DIR, "heatmap_pcos_expression.png"))

# aging heatmap
safe_heatmap(aging_expr, plot_genes, 
             "Aging: Gene Expression by Cell Type (Z-score)", 
             file.path(OUTPUT_DIR, "heatmap_aging_expression.png"))

# cell-type enrichment analysis
# which cell types are most affected?

cat("\nCalculating cell-type enrichment...\n")

# count genes with high expression (>1) per cell type
celltype_enrichment <- combined_expr %>%
  filter(expression > 1) %>%
  group_by(celltype, dataset) %>%
  summarise(
    n_genes_expressed = n_distinct(gene),
    mean_expression = mean(expression),
    .groups = "drop"
  ) %>%
  arrange(desc(n_genes_expressed))

write_csv(
  celltype_enrichment,
  file.path(OUTPUT_DIR, "celltype_enrichment.csv")
)

# barplot of enrichment
p3 <- ggplot(celltype_enrichment, 
             aes(x = reorder(celltype, n_genes_expressed), 
                 y = n_genes_expressed, 
                 fill = dataset)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  theme_minimal() +
  scale_fill_manual(values = c("PCOS" = "#E64B35", "Aging" = "#4DBBD5")) +
  labs(
    title = "Cell-Type Enrichment of Priority Genes",
    x = "Cell Type",
    y = "Number of Genes Expressed (>1)",
    fill = "Dataset"
  )

ggsave(
  file.path(OUTPUT_DIR, "celltype_enrichment_barplot.png"),
  p3,
  width = 8,
  height = 6,
  dpi = 300
)

# identify cell-type markers from priority genes

cat("\nIdentifying cell-type markers...\n")

# define marker as: tau > 0.7 and max expression > 2
# note: tau_pcos might be na if gene wasn't found in pcos, handle that
markers <- tau_scores %>%
  filter(
    (coalesce(tau_PCOS, 0) > 0.7 & coalesce(max_expression_PCOS, 0) > 2) | 
      (coalesce(tau_Aging, 0) > 0.7 & coalesce(max_expression_Aging, 0) > 2)
  ) %>%
  arrange(desc(pmax(coalesce(tau_PCOS,0), coalesce(tau_Aging,0))))

write_csv(
  markers,
  file.path(OUTPUT_DIR, "celltype_specific_markers.csv")
)

# summary report

cat("\n==============================================================================\n")
cat("Cell-Type Specificity Summary\n")
cat("==============================================================================\n\n")

cat("Genes analyzed:", length(top_genes), "\n")
cat("Cell-type specific markers (tau > 0.7):", nrow(markers), "\n\n")

cat("Top 5 Most Enriched Cell Types (PCOS):\n")
print(celltype_enrichment %>% filter(dataset == "PCOS") %>% head(5))

cat("\n\nTop 5 Most Enriched Cell Types (Aging):\n")
print(celltype_enrichment %>% filter(dataset == "Aging") %>% head(5))

cat("\n\nTop 10 Cell-Type Specific Genes:\n")
# select columns safely
cols_to_show <- c("gene", "tau_PCOS", "max_celltype_PCOS", "tau_Aging", "max_celltype_Aging")
cols_present <- intersect(cols_to_show, colnames(markers))
print(markers %>% head(10) %>% select(all_of(cols_present)))

cat("\n\nOutputs saved to:", OUTPUT_DIR, "\n\n")