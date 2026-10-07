#!/usr/bin/env Rscript

# Script to generate a comprehensive transparency table for Tier 1 gene prioritization
# Addressing Major Concern 3 from journal reviewer

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(writexl)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(pheatmap)
})

set.seed(42)

# Define directories
input_dir_prioritization <- "/home/deekshah/mini_project/analysis/09_gene_prioritization/"
input_dir_network <- "/home/deekshah/mini_project/analysis/10_network_analysis/"
input_dir_svg <- "/home/deekshah/mini_project/analysis/16_spatially_variable_genes/"
input_dir_deg <- "/home/deekshah/mini_project/analysis/17_celltype_DEG/"
output_dir <- "/home/deekshah/mini_project/analysis/25_tier1_table/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 1. Load prioritization results
prioritization_file <- file.path(input_dir_prioritization, "Gene_Prioritization_Complete.csv")
if (!file.exists(prioritization_file)) {
  stop(paste("Gene prioritization file not found:", prioritization_file))
}
gene_prior <- read_csv(prioritization_file, show_col_types = FALSE)

# 2. Load PPI network centrality data
network_file <- file.path(input_dir_network, "network_hub_genes.csv")
if (file.exists(network_file)) {
  network_data <- read_csv(network_file, show_col_types = FALSE)
} else {
  warning("Network file not found. Creating empty network data.")
  network_data <- data.frame(gene = character(), degree = numeric(), betweenness = numeric())
}

# 3. Load spatially variable gene data (exclude summary/metadata files)
svg_files <- list.files(input_dir_svg, pattern = "\\.csv$", full.names = TRUE)
svg_files <- svg_files[!grepl("summary|overlap|run_summary|priority", basename(svg_files), ignore.case = TRUE)]
svg_genes <- c()
for (f in svg_files) {
  df <- tryCatch(read_csv(f, show_col_types = FALSE), error = function(e) NULL)
  if (!is.null(df)) {
    # Check common column names for gene symbols
    if ("gene" %in% colnames(df)) {
      svg_genes <- c(svg_genes, df$gene)
    } else if ("X1" %in% colnames(df)) {
      svg_genes <- c(svg_genes, df$X1)
    } else if ("..." %in% colnames(df)) {
      svg_genes <- c(svg_genes, df[[1]])
    } else if (ncol(df) > 0 && is.character(df[[1]])) {
      svg_genes <- c(svg_genes, df[[1]])
    }
  }
}
svg_genes <- unique(svg_genes)

# 4. Create UPDATED prioritization with 6 criteria
# Process network data: normalize degree and betweenness to [0,1], then average
normalize_01 <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return(numeric(0))
  min_x <- min(x, na.rm = TRUE)
  max_x <- max(x, na.rm = TRUE)
  if (max_x == min_x) return(rep(0, length(x)))
  return((x - min_x) / (max_x - min_x))
}

if (nrow(network_data) > 0) {
  network_processed <- network_data %>%
    select(gene, degree, betweenness) %>%
    mutate(
      degree_norm = normalize_01(degree),
      betweenness_norm = normalize_01(betweenness),
      network_score = (degree_norm + betweenness_norm) / 2
    ) %>%
    select(gene, network_score)
} else {
  network_processed <- data.frame(gene = character(), network_score = numeric())
}

# Join and calculate new scores
updated_prior <- gene_prior %>%
  left_join(network_processed, by = "gene") %>%
  mutate(
    network_score = ifelse(is.na(network_score), 0, network_score),
    spatial_score = ifelse(gene %in% svg_genes, 1, 0)
  )

# Ensure concordance_score exists, or derive it if necessary
if (!"concordance_score" %in% names(updated_prior) && "same_direction" %in% names(updated_prior)) {
  updated_prior$concordance_score <- ifelse(updated_prior$same_direction == TRUE, 1, 0)
} else if (!"concordance_score" %in% names(updated_prior)) {
  updated_prior$concordance_score <- 0
}

# Calculate composite score based on defined weights
# Metric 1: Statistical strength (stat_norm) — weight 0.25
# Metric 2: Direction concordance (concordance_score) — weight 0.20
# Metric 3: Pathway centrality (pathway_norm) — weight 0.15
# Metric 4: Effect size magnitude (effect_norm) — weight 0.10
# Metric 5: Network centrality (network_score) — weight 0.20
# Metric 6: Spatial variability (spatial_score) — weight 0.10
updated_prior <- updated_prior %>%
  mutate(
    new_composite_score = 
      (0.25 * stat_norm) + 
      (0.20 * concordance_score) + 
      (0.15 * pathway_norm) + 
      (0.10 * effect_norm) + 
      (0.20 * network_score) + 
      (0.10 * spatial_score)
  )

# 5. Re-assign tiers (Tier 1 = top 20%, Tier 2 = 50-80%, Tier 3 = bottom 50%)
updated_prior <- updated_prior %>%
  arrange(desc(new_composite_score)) %>%
  mutate(
    new_rank = row_number(),
    percentile = 1 - (new_rank - 1) / n(),
    new_tier = case_when(
      percentile > 0.80 ~ "Tier 1",
      percentile > 0.50 ~ "Tier 2",
      TRUE ~ "Tier 3"
    )
  )

# 6. Comprehensive transparency table
transparency_table <- updated_prior %>%
  select(
    gene, 
    new_tier, new_rank, new_composite_score,
    old_tier = tier, old_rank = rank, old_score = final_score,
    metric1_stat_norm = stat_norm, 
    metric2_concordance = concordance_score, 
    metric3_pathway_norm = pathway_norm, 
    metric4_effect_norm = effect_norm, 
    metric5_network_score = network_score, 
    metric6_spatial_score = spatial_score,
    everything()
  )

# 7. Save as CSV and XLSX
write_csv(transparency_table, file.path(output_dir, "Updated_Gene_Prioritization_Transparency.csv"))
write_xlsx(list("Prioritization_Table" = transparency_table), file.path(output_dir, "Updated_Gene_Prioritization_Transparency.xlsx"))

# Generate heatmap of top 30 genes
top30_genes <- transparency_table %>% head(30)
heatmap_data <- top30_genes %>%
  select(gene, metric1_stat_norm, metric2_concordance, metric3_pathway_norm, metric4_effect_norm, metric5_network_score, metric6_spatial_score) %>%
  column_to_rownames("gene")

# TIFF (600 DPI)
tiff(file.path(output_dir, "Top30_Metrics_Heatmap.tiff"), width = 4800, height = 6000, res = 600)
pheatmap(as.matrix(heatmap_data), 
         cluster_cols = FALSE, 
         cluster_rows = TRUE, 
         display_numbers = TRUE,
         main = "Top 30 Prioritized Genes: Transparency Metrics")
dev.off()

# PNG (300 DPI)
png(file.path(output_dir, "Top30_Metrics_Heatmap.png"), width = 2400, height = 3000, res = 300)
pheatmap(as.matrix(heatmap_data), 
         cluster_cols = FALSE, 
         cluster_rows = TRUE, 
         display_numbers = TRUE,
         main = "Top 30 Prioritized Genes: Transparency Metrics")
dev.off()

# 8. Generate a comparison table
comparison_table <- transparency_table %>%
  select(gene, old_tier, new_tier, old_rank, new_rank) %>%
  mutate(
    tier_change = case_when(
      old_tier == new_tier ~ "No Change",
      old_tier == "Tier 1" & new_tier %in% c("Tier 2", "Tier 3") ~ "Downgraded from Tier 1",
      old_tier == "Tier 2" & new_tier == "Tier 3" ~ "Downgraded from Tier 2",
      old_tier %in% c("Tier 2", "Tier 3") & new_tier == "Tier 1" ~ "Upgraded to Tier 1",
      old_tier == "Tier 3" & new_tier == "Tier 2" ~ "Upgraded to Tier 2",
      TRUE ~ "Changed"
    )
  )
write_csv(comparison_table, file.path(output_dir, "Prioritization_Tier_Comparison.csv"))

# 9. Check Runx1 specifically in cell-type DEG files
runx1_report <- list()
if (dir.exists(input_dir_deg)) {
  deg_files <- list.files(input_dir_deg, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  for (f in deg_files) {
    df <- tryCatch(read_csv(f, show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(df) && "gene" %in% colnames(df)) {
      runx1_data <- df %>% filter(gene == "Runx1" | gene == "RUNX1")
      if (nrow(runx1_data) > 0) {
        runx1_report[[basename(f)]] <- runx1_data
      }
    }
  }
}

if (length(runx1_report) > 0) {
  runx1_df <- bind_rows(runx1_report, .id = "source_file")
  write_csv(runx1_df, file.path(output_dir, "Runx1_CellType_DEG_Evidence.csv"))
  cat("Runx1 cell-type DEG evidence found and saved to", file.path(output_dir, "Runx1_CellType_DEG_Evidence.csv"), "\n")
} else {
  cat("No significant Runx1 cell-type DEG evidence found in provided files, or directory empty.\n")
}

cat("Transparency table generation complete. Outputs saved in:", output_dir, "\n")
