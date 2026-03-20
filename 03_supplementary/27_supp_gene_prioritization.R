#!/usr/bin/env Rscript
# script: 09_gene_prioritization_corrected.r
# purpose: multi-metric prioritization of shared pcos-aging genes
# improvements: statistical significance, direction concordance, proper weighting
# reference: crow et al. (2018) nat commun pmid:30209251

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

PROJECT_ROOT <- "E:/Documents/mini_project"

INPUT_DIR <- file.path(
  PROJECT_ROOT,
  "analysis"
)

OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "analysis/09_gene_prioritization"
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# load data

cat("\n==============================================================================\n")
cat("Multi-Metric Gene Prioritization\n")
cat("==============================================================================\n\n")

# load shared degs
shared_deg <- read_csv(
  file.path(INPUT_DIR, "08_venn/Shared_PCOS_Aging_DEGs_full.csv"),
  show_col_types = FALSE
)

# load enrichment results
go_up <- read_csv(
  file.path(INPUT_DIR, "07_go_kegg_overlap/enrichment/GO_BP_shared_up_full.csv"),
  show_col_types = FALSE
)

go_down <- read_csv(
  file.path(INPUT_DIR, "07_go_kegg_overlap/enrichment/GO_BP_shared_down_full.csv"),
  show_col_types = FALSE
)

kegg_up <- read_csv(
  file.path(INPUT_DIR, "07_go_kegg_overlap/enrichment/KEGG_shared_up_full.csv"),
  show_col_types = FALSE
)

kegg_down <- read_csv(
  file.path(INPUT_DIR, "07_go_kegg_overlap/enrichment/KEGG_shared_down_full.csv"),
  show_col_types = FALSE
)

# required object (was missing)
enrich_all <- bind_rows(
  go_up   %>% mutate(source = "GO_up"),
  go_down %>% mutate(source = "GO_down"),
  kegg_up %>% mutate(source = "KEGG_up"),
  kegg_down %>% mutate(source = "KEGG_down")
)

# metric 1: statistical strength score
# combines effect size and significance

cat("Calculating Metric 1: Statistical Strength...\n")

stat_score <- shared_deg %>%
  dplyr::mutate(
    volcano_pcos  = abs(logFC_pcos)  * pmin(-log10(padj_pcos),  10),
    volcano_aging = abs(logFC_aging) * pmin(-log10(padj_aging), 10),
    stat_combined = (volcano_pcos + volcano_aging) / 2,
    stat_norm     = stat_combined / max(stat_combined, na.rm = TRUE)
  ) %>%
  dplyr::select(gene, stat_combined, stat_norm)

# metric 2: direction concordance

cat("Calculating Metric 2: Direction Concordance...\n")

concordance <- shared_deg %>%
  dplyr::mutate(
    same_direction = sign(logFC_pcos) == sign(logFC_aging),
    concordance_score = ifelse(same_direction, 1.0, 0.5),
    direction = dplyr::case_when(
      logFC_pcos  > 0 & logFC_aging > 0 ~ "Both_Up",
      logFC_pcos  < 0 & logFC_aging < 0 ~ "Both_Down",
      TRUE ~ "Discordant"
    )
  ) %>%
  dplyr::select(gene, same_direction, concordance_score, direction)

# metric 3: pathway centrality

cat("Calculating Metric 3: Pathway Centrality...\n")

gene_term_map <- enrich_all %>%
  dplyr::select(Description, geneID, source) %>%
  dplyr::filter(!is.na(geneID), geneID != "") %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(gene = geneID)

pathway_score <- gene_term_map %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    pathway_hits = n(),
    n_go_terms   = sum(grepl("GO", source)),
    n_kegg_terms = sum(grepl("KEGG", source)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    pathway_norm = pathway_hits / max(pathway_hits, na.rm = TRUE)
  )

# metric 4: effect size magnitude

cat("Calculating Metric 4: Effect Size Magnitude...\n")

effect_magnitude <- shared_deg %>%
  dplyr::mutate(
    avg_abs_logFC = (abs(logFC_pcos) + abs(logFC_aging)) / 2,
    effect_norm   = avg_abs_logFC / max(avg_abs_logFC, na.rm = TRUE)
  ) %>%
  dplyr::select(gene, avg_abs_logFC, effect_norm)

# metric 5: expression prevalence
# (columns not present → safe neutral handling)

cat("Calculating Metric 5: Expression Prevalence...\n")

prevalence_score <- shared_deg %>%
  dplyr::mutate(
    expr_prevalence = NA_real_,
    prev_norm = 0
  ) %>%
  dplyr::select(gene, expr_prevalence, prev_norm)

# combine all metrics

cat("\nCombining all metrics...\n")

gene_scores <- shared_deg %>%
  dplyr::select(gene) %>%
  dplyr::left_join(stat_score, by = "gene") %>%
  dplyr::left_join(concordance, by = "gene") %>%
  dplyr::left_join(pathway_score, by = "gene") %>%
  dplyr::left_join(effect_magnitude, by = "gene") %>%
  dplyr::left_join(prevalence_score, by = "gene") %>%
  dplyr::mutate(
    pathway_hits = replace_na(pathway_hits, 0),
    pathway_norm = replace_na(pathway_norm, 0),
    n_go_terms   = replace_na(n_go_terms, 0),
    n_kegg_terms = replace_na(n_kegg_terms, 0)
  )

# final composite score

cat("Calculating final composite scores...\n")

gene_scores <- gene_scores %>%
  dplyr::mutate(
    final_score = (
      0.35 * stat_norm +
        0.25 * concordance_score +
        0.20 * pathway_norm +
        0.15 * effect_norm +
        0.05 * prev_norm
    )
  ) %>%
  dplyr::arrange(desc(final_score))

# tier assignment

cat("Assigning priority tiers...\n")

gene_scores <- gene_scores %>%
  dplyr::mutate(
    tier = dplyr::case_when(
      final_score >= quantile(final_score, 0.80, na.rm = TRUE) ~ "Tier 1: Core Drivers",
      final_score >= quantile(final_score, 0.50, na.rm = TRUE) ~ "Tier 2: Supporting",
      TRUE ~ "Tier 3: Background"
    ),
    rank = dplyr::row_number()
  )

gene_scores <- gene_scores %>%
  dplyr::left_join(
    shared_deg %>% dplyr::select(gene, logFC_pcos, padj_pcos,
                                 logFC_aging, padj_aging),
    by = "gene"
  )

# save outputs

write_csv(
  gene_scores,
  file.path(OUTPUT_DIR, "Gene_Prioritization_Complete.csv")
)

top_genes <- gene_scores %>%
  dplyr::filter(tier == "Tier 1: Core Drivers")

write_csv(
  top_genes,
  file.path(OUTPUT_DIR, "Top_Priority_Genes_Tier1.csv")
)

# visualization

cat("\nGenerating prioritization heatmap...\n")

viz_data <- gene_scores %>%
  dplyr::slice_head(n = 50) %>%
  dplyr::select(gene, stat_norm, concordance_score,
                pathway_norm, effect_norm, prev_norm) %>%
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "metric",
    values_to = "score"
  )

p1 <- ggplot(viz_data, aes(x = metric, y = gene, fill = score)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0.5
  ) +
  theme_minimal()

ggsave(
  file.path(OUTPUT_DIR, "prioritization_heatmap_top50.png"),
  p1,
  width = 8,
  height = 14,
  dpi = 300
)

cat("\nDONE: Gene prioritization completed successfully.\n")