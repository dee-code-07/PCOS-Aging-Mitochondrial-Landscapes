# script: 02_celltype_composition.r
# purpose: cell-type proportion analysis
# output: 600 dpi tiff (nature-ready)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# load data
pcos <- readRDS(
  "E:/Documents/mini_project/scrna/output/07_annotation/pcos/pcos_annotated.rds"
)

aging <- readRDS(
  "E:/Documents/mini_project/scrna/output/07_annotation/aging/aging_annotated.rds"
)

# output directory
outdir <- "E:/Documents/mini_project/analysis/02_celltype_composition"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# nature-style theme
theme_nature <- function() {
  theme_classic(base_size = 12) +
    theme(
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title.x = element_blank(),
      plot.title = element_text(face = "bold", size = 14),
      legend.title = element_blank()
    )
}

# helper function
celltype_composition <- function(seu, condition_label) {
  
  meta <- seu@meta.data %>%
    select(orig.ident, final_celltype)
  
  counts <- meta %>%
    count(orig.ident, final_celltype, name = "n_cells") %>%
    group_by(orig.ident) %>%
    mutate(percent = 100 * n_cells / sum(n_cells)) %>%
    ungroup() %>%
    mutate(condition = condition_label)
  
  return(counts)
}

# run analysis
pcos_comp  <- celltype_composition(pcos,  "PCOS")
aging_comp <- celltype_composition(aging, "Aging")

all_comp <- bind_rows(pcos_comp, aging_comp)

# save tables
write_csv(
  all_comp,
  file.path(outdir, "celltype_composition_percentages.csv")
)

# plot
p <- ggplot(
  all_comp,
  aes(x = orig.ident, y = percent, fill = final_celltype)
) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_wrap(~ condition, scales = "free_x") +
  ylab("Percentage of cells") +
  ggtitle("Cell-type composition in PCOS and Aging ovaries") +
  theme_nature() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# save as 600 dpi tiff
tiff(
  filename = file.path(outdir, "Celltype_composition_PCOS_vs_Aging.tiff"),
  width = 7,
  height = 4,
  units = "in",
  res = 600,
  compression = "lzw"
)
print(p)
dev.off()

message("✅ Cell-type composition analysis completed.")