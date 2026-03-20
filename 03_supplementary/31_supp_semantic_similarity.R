#!/usr/bin/env Rscript
# script: 14_go_semantic_similarity.r
# purpose: cluster genes & apply manual curated names for publication
# method: gosemsim + manual override
# output: functional_gene_modules_curated.csv, similarity_heatmap_curated.png
# date: 2026-02-04

suppressPackageStartupMessages({
  library(GOSemSim)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(dplyr)
  library(readr)
  library(pheatmap)
  library(viridis)
  library(ggplot2)
  library(grid)
  library(RColorBrewer)
  library(stringr)
})

# configuration

PROJECT_ROOT <- "E:/Documents/mini_project"
INPUT_FILE <- file.path(PROJECT_ROOT, "analysis/09_gene_prioritization/Shared_Gene_Prioritization.csv")
OUTPUT_DIR <- file.path(PROJECT_ROOT, "analysis/14_go_semantic")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

TOP_N_GENES <- Inf 
ONTOLOGIES <- c("BP", "MF", "CC") 
SIMILARITY_METHOD <- "Wang"
N_CLUSTERS <- 5

# manual override (publication names)
# based on inspection of your specific gene lists (nr5a2, gja1, ifitm1, etc.)
# these names replace the generic go terms.

CLUSTER_LABELS <- c(
  "1" = "Transcriptional Regulation of Identity",  # Nr5a2, Esr2, Runx1
  "2" = "Immune Response (Interferon Sig.)",       # Ifitm1, Ifitm2
  "3" = "Gap Junctions & Cell Communication",      # Gja1, Dab2
  "4" = "Follicular Remodeling & Activin Reg.",    # Fst, Serpine2
  "5" = "Cytoskeletal Stress & Senescence"         # Tagln, Cnmd
)

# load & map data

cat("\n==============================================================================\n")
cat("GO Semantic Similarity Analysis (Curated Mode)\n")
cat("==============================================================================\n\n")

if(!file.exists(INPUT_FILE)) stop("Input file not found: ", INPUT_FILE)

priority_genes <- read_csv(INPUT_FILE, show_col_types = FALSE)
priority_genes$gene <- trimws(priority_genes$gene)
raw_genes <- priority_genes$gene

# robust id mapping
run_bitr <- function(g, t) tryCatch(bitr(g, fromType=t, toType="ENTREZID", OrgDb=org.Mm.eg.db), error=function(e) NULL)

map_res <- run_bitr(raw_genes, "SYMBOL")
missing <- setdiff(raw_genes, map_res$SYMBOL)
if(length(missing) > 0) {
  res_2 <- run_bitr(str_to_title(missing), "SYMBOL")
  if(!is.null(res_2)) {
    res_2 <- res_2 %>% left_join(data.frame(O=missing, T=str_to_title(missing)), by=c("SYMBOL"="T")) %>% select(SYMBOL=O, ENTREZID)
    map_res <- rbind(map_res, res_2)
  }
}
missing <- setdiff(raw_genes, map_res$SYMBOL)
if(length(missing) > 0) {
  res_3 <- run_bitr(gsub("^Mt-", "mt-", str_to_title(missing)), "SYMBOL")
  if(!is.null(res_3)) {
    res_3 <- res_3 %>% left_join(data.frame(O=missing, M=gsub("^Mt-", "mt-", str_to_title(missing))), by=c("SYMBOL"="M")) %>% select(SYMBOL=O, ENTREZID)
    map_res <- rbind(map_res, res_3)
  }
}

gene_mapping <- map_res[!duplicated(map_res$SYMBOL), ]
entrez_ids <- gene_mapping$ENTREZID
names(entrez_ids) <- gene_mapping$SYMBOL

if(length(entrez_ids) < 5) stop("Too few genes mapped.")

# similarity & clustering

cat("Calculating similarity & clustering...\n")

calculate_go_similarity <- function(genes, ontology, method) {
  all_ids <- as.character(genes)
  full_matrix <- matrix(0, nrow = length(all_ids), ncol = length(all_ids))
  rownames(full_matrix) <- colnames(full_matrix) <- all_ids
  go_data <- godata(OrgDb = org.Mm.eg.db, ont = ontology, computeIC = TRUE)
  sim_result <- mgeneSim(genes = all_ids, semData = go_data, measure = method, combine = "BMA", verbose = FALSE)
  if(is.matrix(sim_result)) {
    valid <- intersect(rownames(sim_result), rownames(full_matrix))
    full_matrix[valid, valid] <- sim_result[valid, valid]
  }
  full_matrix[is.na(full_matrix)] <- 0
  return(full_matrix)
}

sim_list <- list()
for (ont in ONTOLOGIES) {
  tryCatch({ sim_list[[ont]] <- calculate_go_similarity(entrez_ids, ont, SIMILARITY_METHOD) }, error=function(e) cat("Failed:", ont, "\n"))
}

combined_sim <- Reduce("+", sim_list) / length(sim_list)
combined_sim <- (combined_sim + t(combined_sim)) / 2
combined_dist <- as.dist(1 - combined_sim)

hc <- hclust(combined_dist, method = "ward.D2")
raw_clusters <- cutree(hc, k = N_CLUSTERS)

# apply manual names

cat("\nApplying Curated Names...\n")

# create lookup
cluster_df <- data.frame(
  gene_symbol = names(entrez_ids)[match(names(raw_clusters), entrez_ids)],
  entrez_id = names(raw_clusters),
  raw_cluster = raw_clusters,
  stringsAsFactors = FALSE
)

# apply the manual labels (fallback to "cluster x" if label missing)
cluster_df$cluster_name <- sapply(as.character(cluster_df$raw_cluster), function(x) {
  if(x %in% names(CLUSTER_LABELS)) return(CLUSTER_LABELS[[x]])
  return(paste0("Cluster ", x))
})

# merge priority info
cluster_df <- cluster_df %>%
  left_join(priority_genes %>% select(gene, final_score), by = c("gene_symbol" = "gene")) %>%
  arrange(cluster_name, desc(final_score))

write_csv(cluster_df, file.path(OUTPUT_DIR, "functional_gene_modules_curated.csv"))

# visualizations

cat("Generating visualizations...\n")

# 1. heatmap
ann_row <- data.frame(Module = cluster_df$cluster_name)
rownames(ann_row) <- cluster_df$gene_symbol
ids_ordered <- cluster_df$entrez_id
mat_plot <- combined_sim[ids_ordered, ids_ordered]
rownames(mat_plot) <- colnames(mat_plot) <- cluster_df$gene_symbol

unique_mods <- unique(cluster_df$cluster_name)
mod_colors <- setNames(brewer.pal(min(length(unique_mods), 8), "Set2"), unique_mods)
ann_colors <- list(Module = mod_colors)

pheatmap(
  mat_plot, color = viridis(100),
  cluster_rows = FALSE, cluster_cols = FALSE,
  annotation_row = ann_row, annotation_colors = ann_colors,
  show_rownames = ifelse(nrow(mat_plot) > 200, FALSE, TRUE),
  show_colnames = FALSE,
  fontsize_row = ifelse(nrow(mat_plot) > 100, 5, 8),
  main = "Functional Gene Similarity (Curated)",
  filename = file.path(OUTPUT_DIR, "similarity_heatmap_curated.png"),
  width = 16, height = 14
)

# 2. barplot
mod_stats <- cluster_df %>%
  group_by(cluster_name) %>%
  summarise(Avg_Priority = mean(final_score, na.rm=TRUE), Count = n())

p_bar <- ggplot(mod_stats, aes(x=reorder(cluster_name, Avg_Priority), y=Avg_Priority, fill=cluster_name)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=Count), hjust=-0.5) +
  coord_flip() +
  scale_fill_manual(values = mod_colors) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12, face="bold"), legend.position="none") +
  labs(title="Functional Module Importance", subtitle="Based on Multi-Omics Concordance", y="Avg Priority Score", x=NULL)

ggsave(file.path(OUTPUT_DIR, "module_summary_barplot_curated.png"), p_bar, width = 12, height = 7)

cat("\nAnalysis Complete. Curated results saved to:", OUTPUT_DIR, "\n")