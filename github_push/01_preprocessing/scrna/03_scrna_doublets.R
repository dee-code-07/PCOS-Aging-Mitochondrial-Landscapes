#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(scater)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript stage2_doublets_scdblfinder.R <dataset>")
dataset <- args[1]

base_output <- "~/projects/scrna/output"
in_filtered_dir <- file.path(base_output, "02_qc", dataset, "filtered")
out_doublet_dir <- file.path(base_output, "03_doublets", dataset)
dir.create(out_doublet_dir, recursive = TRUE, showWarnings = FALSE)

summary_csv <- file.path(out_doublet_dir, "doublet_summary.csv")

files <- list.files(in_filtered_dir, pattern = "_filtered.rds$", full.names = TRUE)
if (length(files) == 0) stop("No filtered files found.")

summary_df <- data.frame(
  sample = character(),
  before = integer(),
  after = integer(),
  removed = integer(),
  removed_pct = numeric(),
  stringsAsFactors = FALSE
)

for (f in files) {

  sample <- sub("_filtered.rds$", "", basename(f))
  message("---- Running scDblFinder on: ", sample)

  so <- readRDS(f)
  before_n <- ncol(so)

# normalize before converting to sce
  so <- NormalizeData(so, verbose = FALSE)
  DefaultAssay(so) <- "RNA"

  sce <- as.SingleCellExperiment(so)

  sce <- scDblFinder(sce, dbr = 0.04)   # 4% expected doublets

  so$doublet_class <- sce$scDblFinder.class
  so$doublet_score <- sce$scDblFinder.score

# keep only singlets
  so_f <- subset(so, subset = doublet_class == "singlet")
  after_n <- ncol(so_f)

  out_file <- file.path(out_doublet_dir, paste0(sample, "_doubletfiltered.rds"))
  saveRDS(so_f, out_file)

  message("    kept: ", after_n, "/", before_n)

  summary_df[nrow(summary_df)+1,] <- list(
    sample,
    before_n,
    after_n,
    before_n - after_n,
    round((before_n - after_n)/before_n*100, 2)
  )
}

write.table(summary_df, summary_csv, sep=",", row.names=FALSE)
message("DONE. Summary: ", summary_csv)