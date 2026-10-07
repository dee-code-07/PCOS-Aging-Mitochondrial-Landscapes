#!/usr/bin/env Rscript
# script: 10_network_analysis.r
# purpose: protein-protein interaction network analysis of prioritized genes
# databases: string v12.0
# output: network graphs, hub gene identification, sub-modules
# reference: szklarczyk et al. (2021) nucleic acids res pmid:33237311

suppressPackageStartupMessages({
  library(STRINGdb)
  library(igraph)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggraph)
  library(tidygraph)
})

# configuration

PROJECT_ROOT <- "E:/Documents/mini_project"

INPUT_FILE <- file.path(
  PROJECT_ROOT,
  "analysis/09_gene_prioritization/Gene_Prioritization_Complete.csv"
)

OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "analysis/10_network_analysis"
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# string parameters
SPECIES_ID <- 10090  # Mus musculus (mouse)
SCORE_THRESHOLD <- 400  # Medium confidence (0-1000 scale)

# network parameters
MIN_CLUSTER_SIZE <- 5
TOP_N_GENES <- 100  # Analyze top 100 prioritized genes

# load data

cat("\n==============================================================================\n")
cat("STRING Protein-Protein Interaction Network Analysis\n")
cat("==============================================================================\n\n")

gene_priority <- read_csv(INPUT_FILE, show_col_types = FALSE)

# select top genes
top_genes <- gene_priority %>%
  dplyr::arrange(desc(final_score)) %>%
  head(TOP_N_GENES) %>%
  dplyr::pull(gene)

cat("Analyzing top", TOP_N_GENES, "genes\n")
cat("Score threshold:", SCORE_THRESHOLD, "(medium confidence)\n\n")

# connect to string database

cat("Connecting to STRING database...\n")
string_db <- STRINGdb$new(
  version = "12.0",
  species = SPECIES_ID,
  score_threshold = SCORE_THRESHOLD,
  network_type = "full",
  input_directory = ""
)

# map genes to string ids

cat("Mapping genes to STRING identifiers...\n")

gene_df <- data.frame(
  gene = top_genes,
  stringsAsFactors = FALSE
)

mapped <- string_db$map(
  gene_df,
  "gene",
  removeUnmappedRows = TRUE
)

cat(sprintf(
  "  Mapped: %d/%d genes (%.1f%%)\n",
  nrow(mapped),
  length(top_genes),
  100 * nrow(mapped) / length(top_genes)
))

write_csv(
  mapped,
  file.path(OUTPUT_DIR, "string_gene_mapping.csv")
)

# retrieve interactions

cat("\nRetrieving protein-protein interactions...\n")

interactions <- string_db$get_interactions(mapped$STRING_id)

cat(sprintf("  Found %d interactions\n", nrow(interactions)))

interactions_genes <- interactions %>%
  dplyr::left_join(
    mapped %>% dplyr::select(STRING_id, from_gene = gene),
    by = c("from" = "STRING_id")
  ) %>%
  dplyr::left_join(
    mapped %>% dplyr::select(STRING_id, to_gene = gene),
    by = c("to" = "STRING_id")
  ) %>%
  dplyr::filter(!is.na(from_gene), !is.na(to_gene)) %>%
  dplyr::select(from_gene, to_gene, combined_score, everything())

write_csv(
  interactions_genes,
  file.path(OUTPUT_DIR, "ppi_interactions.csv")
)

# create igraph network object

cat("\nBuilding network graph...\n")

edges <- interactions_genes %>%
  dplyr::select(from = from_gene, to = to_gene, weight = combined_score)

g <- graph_from_data_frame(
  d = edges,
  directed = FALSE,
  vertices = unique(c(edges$from, edges$to))
)

node_attrs <- gene_priority %>%
  dplyr::filter(gene %in% V(g)$name) %>%
  dplyr::select(gene, final_score, tier)

V(g)$priority_score <- node_attrs$final_score[match(V(g)$name, node_attrs$gene)]
V(g)$tier <- node_attrs$tier[match(V(g)$name, node_attrs$gene)]

V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g)
V(g)$closeness <- closeness(g)
V(g)$eigen_centrality <- eigen_centrality(g)$vector

# identify hub genes

cat("\nIdentifying hub genes...\n")

hub_genes <- data.frame(
  gene = V(g)$name,
  degree = V(g)$degree,
  betweenness = V(g)$betweenness,
  closeness = V(g)$closeness,
  eigen_centrality = V(g)$eigen_centrality,
  priority_score = V(g)$priority_score,
  tier = V(g)$tier,
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(desc(degree))

hub_threshold <- quantile(hub_genes$degree, 0.90, na.rm = TRUE)

hub_genes <- hub_genes %>%
  dplyr::mutate(is_hub = degree >= hub_threshold)

write_csv(
  hub_genes,
  file.path(OUTPUT_DIR, "network_hub_genes.csv")
)

cat(sprintf(
  "  Identified %d hub genes (top 10%% by degree)\n",
  sum(hub_genes$is_hub)
))

# community detection (modules)

cat("\nDetecting network modules...\n")

communities <- cluster_louvain(g)
V(g)$module <- membership(communities)

module_summary <- data.frame(
  gene = V(g)$name,
  module = V(g)$module,
  priority_score = V(g)$priority_score,
  tier = V(g)$tier
) %>%
  dplyr::group_by(module) %>%
  dplyr::summarise(
    n_genes = n(),
    mean_priority = mean(priority_score, na.rm = TRUE),
    genes = paste(gene, collapse = ", "),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_genes >= MIN_CLUSTER_SIZE) %>%
  dplyr::arrange(desc(mean_priority))

write_csv(
  module_summary,
  file.path(OUTPUT_DIR, "network_modules.csv")
)

# visualization

cat("\nGenerating network visualizations...\n")

tg <- as_tbl_graph(g)

p1 <- ggraph(tg, layout = "fr") +
  geom_edge_link(aes(alpha = weight), color = "gray70", show.legend = FALSE) +
  geom_node_point(aes(size = degree, color = priority_score), alpha = 0.8) +
  geom_node_text(
    aes(label = ifelse(degree > hub_threshold, name, "")),
    repel = TRUE, size = 3, fontface = "bold"
  ) +
  scale_color_gradient(low = "lightblue", high = "red") +
  scale_size_continuous(range = c(2, 10)) +
  theme_graph()

ggsave(
  file.path(OUTPUT_DIR, "ppi_network_full.png"),
  p1,
  width = 14,
  height = 12,
  dpi = 300
)

# export for cytoscape

cytoscape_nodes <- data.frame(
  id = V(g)$name,
  gene = V(g)$name,
  degree = V(g)$degree,
  betweenness = V(g)$betweenness,
  priority_score = V(g)$priority_score,
  tier = V(g)$tier,
  module = V(g)$module,
  is_hub = V(g)$name %in% hub_genes$gene[hub_genes$is_hub]
)

write_csv(
  cytoscape_nodes,
  file.path(OUTPUT_DIR, "cytoscape_nodes.csv")
)

cytoscape_edges <- interactions_genes %>%
  dplyr::select(
    source = from_gene,
    target = to_gene,
    weight = combined_score
  )

write_csv(
  cytoscape_edges,
  file.path(OUTPUT_DIR, "cytoscape_edges.csv")
)

cat("\nDONE: Network analysis completed successfully.\n")