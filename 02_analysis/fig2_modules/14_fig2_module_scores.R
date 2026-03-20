# script: 03_module_scores.r
# purpose: add biological module scores
# output: seurat objects with module scores

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
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
outdir <- "E:/Documents/mini_project/analysis/03_modules"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# module gene sets
module_dir <- "E:/Documents/mini_project/scrna/resources/modules"

modules <- list(
  Mitochondrial = "mitochondrial_module.csv",
  Senescence    = "senescence_SASP_module.csv",
  Oxidative     = "oxidative_stress_module.csv",
  Inflammation  = "inflammation_module.csv"
)

# scoring function
score_modules <- function(seu) {
  
  genes_present <- toupper(rownames(seu))
  
  for (m in names(modules)) {
    
    genes <- read_csv(
      file.path(module_dir, modules[[m]]),
      show_col_types = FALSE
    )[[1]] |>
      toupper() |>
      intersect(genes_present)
    
    if (length(genes) < 5) next
    
    seu <- AddModuleScore(
      seu,
      features = list(genes),
      name = paste0(m, "_Score"),
      assay = "SCT",
      search = FALSE
    )
  }
  
  return(seu)
}

# run module scoring
pcos  <- score_modules(pcos)
aging <- score_modules(aging)

# save outputs
saveRDS(
  pcos,
  file.path(outdir, "pcos_with_modules.rds")
)

saveRDS(
  aging,
  file.path(outdir, "aging_with_modules.rds")
)

message("✅ Module scoring completed successfully.")