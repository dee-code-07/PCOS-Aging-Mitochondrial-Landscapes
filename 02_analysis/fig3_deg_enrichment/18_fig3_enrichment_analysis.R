# 02_run_go_kegg_shared.r
# go bp + kegg enrichment with full annotation

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(readr)
})

BASE_DIR <- "E:/Documents/mini_project/analysis/07_go_kegg_overlap"
INPUT_DIR <- file.path(BASE_DIR, "input")
OUT_DIR   <- file.path(BASE_DIR, "enrichment")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# load shared deg tables
shared_up <- read_csv(
  file.path(INPUT_DIR, "Shared_DEGs_up_full.csv"),
  show_col_types = FALSE
)

shared_down <- read_csv(
  file.path(INPUT_DIR, "Shared_DEGs_down_full.csv"),
  show_col_types = FALSE
)

# map gene symbols → entrez ids
map_genes <- function(df) {
  mapped <- bitr(
    df$gene,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Mm.eg.db
  )
  inner_join(df, mapped, by = c("gene" = "SYMBOL"))
}

shared_up_mapped   <- map_genes(shared_up)
shared_down_mapped <- map_genes(shared_down)

# go bp enrichment
run_go <- function(entrez_ids) {
  enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
}

ego_up   <- run_go(shared_up_mapped$ENTREZID)
ego_down <- run_go(shared_down_mapped$ENTREZID)

# go cc enrichment
run_go_cc <- function(entrez_ids) {
  enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Mm.eg.db,
    keyType       = "ENTREZID",
    ont           = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE
  )
}

ego_cc_up   <- run_go_cc(shared_up_mapped$ENTREZID)
ego_cc_down <- run_go_cc(shared_down_mapped$ENTREZID)

# kegg enrichment
run_kegg <- function(entrez_ids) {
  enrichKEGG(
    gene          = entrez_ids,
    organism      = "mmu",
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1
  )
}

ekegg_up   <- run_kegg(shared_up_mapped$ENTREZID)
ekegg_down <- run_kegg(shared_down_mapped$ENTREZID)

# save full enrichment tables
write_csv(as.data.frame(ego_up),
          file.path(OUT_DIR, "GO_BP_shared_up_full.csv"))

write_csv(as.data.frame(ego_down),
          file.path(OUT_DIR, "GO_BP_shared_down_full.csv"))

write_csv(as.data.frame(ego_cc_up),
          file.path(OUT_DIR, "GO_CC_shared_up_full.csv"))

write_csv(as.data.frame(ego_cc_down),
          file.path(OUT_DIR, "GO_CC_shared_down_full.csv"))


write_csv(as.data.frame(ekegg_up),
          file.path(OUT_DIR, "KEGG_shared_up_full.csv"))

write_csv(as.data.frame(ekegg_down),
          file.path(OUT_DIR, "KEGG_shared_down_full.csv"))

message("✅ GO + KEGG enrichment completed.")
message("Files written to: ", OUT_DIR)