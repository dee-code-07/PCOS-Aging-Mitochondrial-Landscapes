# script: 01_umap_overview.r
# purpose: umap visualization for pcos and aging scrna-seq

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# load data
pcos <- readRDS(
  "E:/Documents/mini_project/scrna/output/07_annotation/pcos/pcos_annotated.rds"
)

aging <- readRDS(
  "E:/Documents/mini_project/scrna/output/07_annotation/aging/aging_annotated.rds"
)

DefaultAssay(pcos)  <- "SCT"
DefaultAssay(aging) <- "SCT"

# output directory
outdir <- "E:/Documents/mini_project/analysis/01_umap_visualization"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# nature-style theme (minimal)
theme_nature <- function() {
  theme_classic(base_size = 12) +
    theme(
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_blank(),
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", size = 14)
    )
}

# umap: clusters
p1 <- DimPlot(
  pcos,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) +
  ggtitle("PCOS – Clusters") +
  theme_nature()

p2 <- DimPlot(
  aging,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
) +
  ggtitle("Aging – Clusters") +
  theme_nature()

ggsave(
  filename = file.path(outdir, "UMAP_clusters_PCOS_vs_Aging.tiff"),
  plot = p1 + p2,
  width = 12,
  height = 6,
  dpi = 600
)

# umap: cell types
p3 <- DimPlot(
  pcos,
  reduction = "umap",
  group.by = "final_celltype",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("PCOS – Cell Types") +
  theme_nature()

p4 <- DimPlot(
  aging,
  reduction = "umap",
  group.by = "final_celltype",
  label = TRUE,
  repel = TRUE
) +
  ggtitle("Aging – Cell Types") +
  theme_nature()

ggsave(
  filename = file.path(outdir, "UMAP_celltypes_PCOS_vs_Aging.tiff"),
  plot = p3 + p4,
  width = 14,
  height = 6,
  dpi = 600
)

message("✅ UMAP overview completed successfully.")