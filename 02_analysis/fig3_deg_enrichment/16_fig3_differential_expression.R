# script: 06_pcos_differential_expression.r
# purpose: pcos case vs control de (sct-safe)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
})

# load pcos object with modules
pcos <- readRDS(
  "E:/Documents/mini_project/analysis/03_modules/pcos_with_modules.rds"
)

DefaultAssay(pcos) <- "SCT"

# define condition using sample_id
pcos$condition <- ifelse(
  grepl("case", pcos$sample_id, ignore.case = TRUE),
  "PCOS_case",
  "PCOS_control"
)

# sanity check
table(pcos$condition)

# prepare sct object for de  <<< this is the fix
pcos <- PrepSCTFindMarkers(pcos)

# set identities
Idents(pcos) <- "condition"

# differential expression
deg_pcos <- FindMarkers(
  pcos,
  ident.1 = "PCOS_case",
  ident.2 = "PCOS_control",
  test.use = "wilcox",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# save results
outdir <- "E:/Documents/mini_project/scrna/output/08_de/pcos"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

deg_pcos <- deg_pcos %>%
  rownames_to_column("gene")

write_csv(
  deg_pcos,
  file.path(outdir, "DE_pcos_case_vs_ctrl.csv")
)

message("✅ PCOS differential expression completed successfully.")