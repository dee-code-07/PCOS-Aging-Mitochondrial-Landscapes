#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript stage5_find_markers.R <dataset>")
dataset <- args[1]

# paths
base_out <- "~/projects/scrna/output"
integrate_dir <- file.path(base_out, "05_integration", dataset)
annot_out_dir  <- file.path(base_out, "06_markers", dataset)
plots_dir      <- file.path(annot_out_dir, "plots")

dir.create(annot_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# parameters (tweak if needed)
params <- list(
  dataset = dataset,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  top_n = 10,
  marker_timestamp = as.character(Sys.time())
)

# load integrated object
integrated_file <- file.path(integrate_dir, paste0(dataset, "_integrated.rds"))
if (!file.exists(integrated_file)) stop("integrated RDS not found: ", integrated_file)
message("loading integrated object: ", integrated_file)
so <- readRDS(integrated_file)

# ensure rna layers collapsed for safety earlier, but markers use sct
DefaultAssay(so) <- params$assay
so <- PrepSCTFindMarkers(so)

message("running FindAllMarkers (this can be slow)...")
markers <- FindAllMarkers(so,
                          assay = params$assay,
                          only.pos = params$only.pos,
                          min.pct = params$min.pct,
                          logfc.threshold = params$logfc.threshold)

# write outputs
write_csv(markers, file.path(annot_out_dir, "markers_all.csv"))
topn <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = params$top_n)
write_csv(topn, file.path(annot_out_dir, "markers_top10.csv"))

# umap by cluster
p1 <- DimPlot(so, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle(paste(dataset, "clusters"))
ggsave(file.path(plots_dir, "umap_clusters.pdf"), p1, width = 8, height = 6)

# top markers heatmap (top 5 per cluster; will skip if huge)
top5 <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 5)
genes_heat <- unique(top5$gene)
if (length(genes_heat) > 1 && length(genes_heat) <= 500) {
  pheat <- DoHeatmap(so, features = genes_heat, assay = params$assay, label = FALSE) + ggtitle("top5 markers heatmap")
  ggsave(file.path(plots_dir, "markers_top5_heatmap.pdf"), pheat, width = 10, height = 8)
}

# save parameter log and small summary
params_df <- tibble::enframe(params, name = "param", value = "value")
write_csv(params_df, file.path(annot_out_dir, "stage5_parameters.csv"))

summary_df <- tibble(
  dataset = dataset,
  n_clusters = length(unique(markers$cluster)),
  total_markers = nrow(markers),
  timestamp = as.character(Sys.time())
)
write_csv(summary_df, file.path(annot_out_dir, "stage5_summary.csv"))

# save object with marker flag
so$markers_called <- TRUE
saveRDS(so, file.path(annot_out_dir, paste0(dataset, "_markers_called.rds")))

message("stage5 complete -> ", annot_out_dir)
