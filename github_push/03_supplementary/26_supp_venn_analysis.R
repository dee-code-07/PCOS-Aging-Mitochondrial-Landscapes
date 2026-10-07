# script: 08_pcos_aging_venn.r
# purpose: publication-grade pcos vs aging venn (nature style)

suppressPackageStartupMessages({
  library(VennDiagram)
  library(grid)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
})

# load degs
deg_pcos <- readr::read_csv(
  "E:/Documents/mini_project/scrna/output/08_de/pcos/DE_pcos_case_vs_ctrl.csv",
  show_col_types = FALSE
)

deg_aging <- readr::read_csv(
  "E:/Documents/mini_project/scrna/output/08_de/aging/DE_aging_aged_vs_young.csv",
  show_col_types = FALSE
)

# thresholds
padj_cutoff <- 0.05
lfc_cutoff  <- 0.25

# significant genes (pcos)
pcos_sig <- deg_pcos %>%
  dplyr::filter(
    p_val_adj < padj_cutoff,
    abs(avg_log2FC) > lfc_cutoff
  ) %>%
  dplyr::select(
    gene,
    logFC_pcos = avg_log2FC,
    padj_pcos  = p_val_adj
  )

# significant genes (aging)
aging_sig <- deg_aging %>%
  dplyr::filter(
    p_val_adj < padj_cutoff,
    abs(avg_log2FC) > lfc_cutoff
  ) %>%
  dplyr::select(
    gene,
    logFC_aging = avg_log2FC,
    padj_aging  = p_val_adj
  )

# shared genes with direction consistency
shared_df <- dplyr::inner_join(
  pcos_sig,
  aging_sig,
  by = "gene"
) %>%
  dplyr::mutate(
    dir_pcos  = ifelse(logFC_pcos  > 0, "Up", "Down"),
    dir_aging = ifelse(logFC_aging > 0, "Up", "Down"),
    same_direction = dir_pcos == dir_aging,
    combined_effect = (abs(logFC_pcos) + abs(logFC_aging)) / 2
  ) %>%
  dplyr::filter(same_direction)

shared_genes <- shared_df$gene

# top shared genes (combined effect size)
top_shared <- shared_df %>%
  dplyr::arrange(dplyr::desc(combined_effect)) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(gene)

# format gene text block
formatted_gene_text <- stringr::str_wrap(
  paste(top_shared, collapse = ", "),
  width = 60
)

# output setup
outdir <- "E:/Documents/mini_project/analysis/08_venn"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tiff(
  filename = file.path(outdir, "Venn_PCOS_vs_Aging_Nature.tiff"),
  width = 6,
  height = 6,
  units = "in",
  res = 600,
  compression = "lzw"
)

# draw venn (style unchanged)
grid.newpage()

venn.plot <- VennDiagram::draw.pairwise.venn(
  area1 = nrow(pcos_sig),
  area2 = nrow(aging_sig),
  cross.area = length(shared_genes),
  category = c("PCOS", "Aging"),
  
  fill = c("#D55E00", "#0072B2"),
  alpha = c(0.6, 0.6),
  lty = "blank",
  
  fontfamily = "sans",
  cat.fontfamily = "sans",
  
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  cat.cex = 1.2,
  cat.fontface = "bold",
  cex = 1.5,
  
  margin = 0.10,
  scaled = FALSE
)

grid::grid.draw(venn.plot)

# leader lines
grid::grid.lines(
  x = c(0.50, 0.20),
  y = c(0.42, 0.18),
  gp = grid::gpar(col = "black", lwd = 1)
)

grid::grid.lines(
  x = c(0.50, 0.80),
  y = c(0.42, 0.18),
  gp = grid::gpar(col = "black", lwd = 1)
)

# shared genes text block
grid::grid.text(
  label = formatted_gene_text,
  x = 0.50,
  y = 0.12,
  just = "center",
  gp = grid::gpar(
    fontsize = 9,
    fontface = "bold",
    fontfamily = "sans",
    col = "#404040",
    lineheight = 1.2
  )
)

# finalize
dev.off()

# save full shared deg table (important for supplementary)
readr::write_csv(
  shared_df,
  file.path(outdir, "Shared_PCOS_Aging_DEGs_full.csv")
)

message("✅ Nature-style Venn generated with direction-consistent top shared genes.")