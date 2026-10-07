#!/usr/bin/env Rscript
# script: 12_tf_enrichment.r
# project: mitochondrial dysfunction in pcos and ovarian aging
# purpose: identify transcription factors regulating priority genes
# databases: consensus (human/mouse) for maximum coverage
# method: enrichr api
# date: 2026-02-04

suppressPackageStartupMessages({
  library(enrichR)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)    # Required for unnest()
  library(stringr)  # Required for str_trunc()
  library(grid)     # Required for margins
})

# configuration & paths

PROJECT_ROOT <- "E:/Documents/mini_project"
INPUT_FILE <- file.path(PROJECT_ROOT, "analysis/09_gene_prioritization/Top_Priority_Genes_Tier1.csv")
OUTPUT_DIR <- file.path(PROJECT_ROOT, "analysis/12_tf_enrichment")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# enrichr tf databases (consensus human/mouse for robustness)
TF_DATABASES <- c(
  "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
  "TRRUST_Transcription_Factors_2019",
  "ChEA_2022",
  "TRANSFAC_and_JASPAR_PWMs"
)

TOP_N_GENES <- 100  # Analyze top 100 genes
TOP_N_TFS <- 30     # Report top 30 TFs

# load data

cat("\n==============================================================================\n")
cat("Transcription Factor Enrichment Analysis\n")
cat("==============================================================================\n\n")

if(!file.exists(INPUT_FILE)) stop("Input file not found: ", INPUT_FILE)

priority_genes <- read_csv(INPUT_FILE, show_col_types = FALSE)
# ensure we don't try to take more genes than exist
n_genes <- min(TOP_N_GENES, nrow(priority_genes))
gene_list <- priority_genes$gene[1:n_genes]

cat("Analyzing", length(gene_list), "prioritized genes\n")
cat("Using databases:", paste(TF_DATABASES, collapse = ", "), "\n\n")

# set enrichr organism

# check internet connection
if(!curl::has_internet()) stop("No internet connection! Enrichr requires active internet.")

setEnrichrSite("Enrichr")  # Default: human/mouse combined

# run enrichment

cat("\nRunning TF enrichment...\n")

enrichment_results <- tryCatch({
  enrichr(genes = gene_list, databases = TF_DATABASES)
}, error = function(e) stop("Enrichr API failed: ", e$message))

# process results

# combine all results
all_results <- bind_rows(
  lapply(names(enrichment_results), function(db) {
    df <- enrichment_results[[db]]
    if(!is.null(df) && nrow(df) > 0) {
      df %>% mutate(database = db)
    } else {
      NULL
    }
  })
)

if(nrow(all_results) == 0) stop("No significant enrichment results found.")

# filter significant results
sig_results <- all_results %>%
  filter(Adjusted.P.value < 0.05) %>%
  arrange(Adjusted.P.value)

write_csv(sig_results, file.path(OUTPUT_DIR, "tf_enrichment_significant.csv"))

# top tfs
top_tfs <- sig_results %>%
# clean names for grouping (remove parenthesis content and species tags)
  mutate(TF_Clean = gsub("\\s*\\(.*$", "", Term)) %>%
  mutate(TF_Clean = gsub("\\s+human|\\s+mouse", "", TF_Clean, ignore.case=TRUE)) %>%
  group_by(TF_Clean) %>%
  summarise(
    Term_Full = first(Term),
    n_databases = n_distinct(database),
    min_pvalue = min(Adjusted.P.value),
    databases = paste(unique(database), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(desc(n_databases), min_pvalue) %>%
  head(TOP_N_TFS)

write_csv(top_tfs, file.path(OUTPUT_DIR, "top_transcription_factors.csv"))

# visualization 1: tf enrichment barplot

cat("Generating visualizations...\n")

# extract tf name and format for plotting
plot_data <- sig_results %>%
  head(TOP_N_TFS) %>%
  mutate(
    TF = gsub("\\s*\\(.*$", "", Term),  # Remove parentheses part
    TF = gsub("\\s+human|\\s+mouse", "", TF, ignore.case=TRUE), # Remove species tags
    TF = str_trunc(TF, 30), # Truncate long names to prevent cutoff
    neg_log10_p = -log10(Adjusted.P.value)
  )

p1 <- ggplot(plot_data, aes(x = reorder(TF, neg_log10_p), y = neg_log10_p, fill = database)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = paste0("Top ", TOP_N_TFS, " Enriched Transcription Factors"),
    x = NULL,
    y = "-log10(Adjusted P-value)",
    fill = "Database"
  ) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(t=10, r=10, b=10, l=10, unit="pt") # Extra margin for text
  )

ggsave(
  file.path(OUTPUT_DIR, "tf_enrichment_barplot.png"),
  p1,
  width = 12, # Increased width to fit text
  height = 8,
  dpi = 300
)

# visualization 2: dot plot (multi-database)

# create matrix: tf x database
tf_matrix <- sig_results %>%
  head(50) %>%
  mutate(
    TF = gsub("\\s*\\(.*$", "", Term),
    TF = gsub("\\s+human|\\s+mouse", "", TF, ignore.case=TRUE),
    neg_log10_p = -log10(Adjusted.P.value)
  ) %>%
  select(TF, database, neg_log10_p)

p2 <- ggplot(tf_matrix, aes(x = database, y = TF, size = neg_log10_p, color = neg_log10_p)) +
  geom_point() +
  scale_size_continuous(range = c(2, 10), name = "-log10(P)") +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(P)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90"),
    plot.margin = margin(t=10, r=10, b=10, l=10, unit="pt")
  ) +
  labs(
    title = "TF Enrichment Across Databases",
    x = NULL,
    y = NULL
  )

ggsave(
  file.path(OUTPUT_DIR, "tf_enrichment_dotplot.png"),
  p2,
  width = 11,
  height = 12,
  dpi = 300
)

# extract target genes for top tfs

cat("\nExtracting target genes for top TFs...\n")

# get genes for each top tf
tf_targets <- sig_results %>%
  filter(Term %in% top_tfs$Term_Full) %>%
  select(Term, Genes, database, Adjusted.P.value) %>%
  mutate(
    TF = gsub("\\s*\\(.*$", "", Term),
    target_genes = strsplit(Genes, ";")
  ) %>%
  unnest(target_genes) %>%
  select(TF, target_gene = target_genes, database, p_value = Adjusted.P.value)

write_csv(tf_targets, file.path(OUTPUT_DIR, "tf_target_genes.csv"))

# identify key regulatory tfs

# define key tfs based on biological relevance (pcos/aging)
# reference: plu et al. (2020) front endocrinol
key_tfs_reference <- c(
  "FOXO1", "FOXO3",     # Aging, apoptosis
  "NFKB1", "RELA",      # NF-kB inflammation
  "NFE2L2",             # NRF2 oxidative stress
  "TP53",               # p53 senescence
  "STAT3",              # JAK-STAT inflammation
  "ESR1", "ESR2",       # Estrogen receptor
  "AR",                 # Androgen receptor (PCOS)
  "PPARGC1A",           # PGC-1a mitochondria
  "TFAM",               # Mitochondrial biogenesis
  "SREBF1", "SREBF2"    # Lipid metabolism
)

cat("Checking for key regulatory TFs...\n")

# check which key tfs are enriched
key_tf_enriched <- sig_results %>%
  mutate(TF_Clean = toupper(gsub("\\s*\\(.*$", "", Term))) %>% # Clean name to match reference
  filter(TF_Clean %in% toupper(key_tfs_reference))

write_csv(key_tf_enriched, file.path(OUTPUT_DIR, "key_regulatory_tfs.csv"))

# summary
cat("\n==============================================================================\n")
cat("Transcription Factor Enrichment Summary\n")
cat("==============================================================================\n\n")
cat("Total TFs identified:", n_distinct(sig_results$Term), "\n")
cat("Significant (FDR < 0.05):", nrow(sig_results), "\n")
cat("Key disease-relevant TFs found:", nrow(key_tf_enriched), "\n\n")

cat("Top 10 TFs:\n")
print(top_tfs %>% select(TF_Clean, n_databases, min_pvalue) %>% head(10))

cat("\nKey Regulatory TFs (Disease-relevant):\n")
if (nrow(key_tf_enriched) > 0) {
  print(key_tf_enriched %>% select(Term, Adjusted.P.value, database))
} else {
  cat("  None found in reference list\n")
}

cat("\nOutputs saved to:", OUTPUT_DIR, "\n\n")