# 03_plot_go_kegg_overlap_simple.r
# simple go + kegg dotplots and barplots (600 dpi)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(readr)
})

BASE_DIR <- "E:/Documents/mini_project/analysis/07_go_kegg_overlap"
IN_DIR   <- file.path(BASE_DIR, "enrichment")
OUT_DIR  <- file.path(BASE_DIR, "figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# helper: csv → enrichresult
csv_to_enrich <- function(df, ontology) {
  
  if (!"ID" %in% colnames(df)) {
    df$ID <- paste0("term_", seq_len(nrow(df)))
  }
  
  new(
    "enrichResult",
    result        = df,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    qvalueCutoff  = 1,
    ontology      = ontology,
    organism      = "mouse"
  )
}

# load enrichment csvs
go_up_df    <- read_csv(file.path(IN_DIR, "GO_BP_shared_up_full.csv"),
                        show_col_types = FALSE)
go_down_df  <- read_csv(file.path(IN_DIR, "GO_BP_shared_down_full.csv"),
                        show_col_types = FALSE)
kegg_up_df  <- read_csv(file.path(IN_DIR, "KEGG_shared_up_full.csv"),
                        show_col_types = FALSE)
kegg_down_df<- read_csv(file.path(IN_DIR, "KEGG_shared_down_full.csv"),
                        show_col_types = FALSE)

ego_up    <- csv_to_enrich(go_up_df,   "BP")
ego_down  <- csv_to_enrich(go_down_df, "BP")
ekegg_up  <- csv_to_enrich(kegg_up_df, "KEGG")
ekegg_down<- csv_to_enrich(kegg_down_df,"KEGG")

# tiff saver (600 dpi)
options(bitmapType = "cairo")

save_tiff <- function(p, fname, w = 8, h = 10) {
  tiff(file.path(OUT_DIR, fname),
       width = w, height = h, units = "in",
       res = 600, compression = "lzw", type = "cairo")
  print(p)
  dev.off()
  message("Saved: ", fname)
}

# go bp — up
save_tiff(
  dotplot(ego_up, showCategory = 20),
  "GO_BP_shared_up_dotplot_600dpi.tiff"
)

save_tiff(
  barplot(ego_up, showCategory = 20),
  "GO_BP_shared_up_barplot_600dpi.tiff"
)

# go bp — down
save_tiff(
  dotplot(ego_down, showCategory = 20),
  "GO_BP_shared_down_dotplot_600dpi.tiff"
)

save_tiff(
  barplot(ego_down, showCategory = 20),
  "GO_BP_shared_down_barplot_600dpi.tiff"
)

# kegg — up
save_tiff(
  dotplot(ekegg_up, showCategory = 20),
  "KEGG_shared_up_dotplot_600dpi.tiff"
)

save_tiff(
  barplot(ekegg_up, showCategory = 20),
  "KEGG_shared_up_barplot_600dpi.tiff"
)

# kegg — down
save_tiff(
  dotplot(ekegg_down, showCategory = 20),
  "KEGG_shared_down_dotplot_600dpi.tiff"
)

save_tiff(
  barplot(ekegg_down, showCategory = 20),
  "KEGG_shared_down_barplot_600dpi.tiff"
)

message("✅ Simple GO + KEGG plots generated (dotplots + barplots).")