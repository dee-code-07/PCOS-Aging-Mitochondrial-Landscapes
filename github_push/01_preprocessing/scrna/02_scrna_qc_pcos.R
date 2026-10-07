#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(readr); library(tibble)
})

dataset <- "pcos"
base_parent <- "~/projects/scrna/data/GSE268919_PCOS_SEQ"

out_base <- "~/projects/scrna/output"
out_load <- file.path(out_base, "01_load", dataset)
out_qc <- file.path(out_base, "02_qc", dataset)
out_qc_plots <- file.path(out_qc, "plots")
out_filtered <- file.path(out_qc, "filtered")
dirs <- c(out_load, out_qc, out_qc_plots, out_filtered)
sapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

params <- list(
  dataset = dataset,
  min_features = 300,
  min_counts = 1000,
  percent_mt_max = 20,
  qc_timestamp = as.character(Sys.time())
)

h5_files <- list.files(base_parent, pattern = "\\.h5$", full.names = TRUE)
samples <- c()

for (f in h5_files) {
  samp <- sub("\\.h5$", "", basename(f)); samples <- c(samples, samp)
  message("loading h5: ", samp)
  raw <- Read10X_h5(f)
  if (is.list(raw)) {
    if ("Gene Expression" %in% names(raw)) counts <- raw[["Gene Expression"]]
    else if ("RNA" %in% names(raw)) counts <- raw[["RNA"]]
    else stop("no rna/gene expression matrix in: ", f)
  } else counts <- raw

  so <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200, project = "PCOS")
  so$sample_id <- samp
  so$percent.mt <- PercentageFeatureSet(so, pattern = "^mt-|^MT-|^mt")
  saveRDS(so, file.path(out_load, paste0(samp, "_raw.rds")))
  message("saved raw rds: ", file.path(out_load, paste0(samp, "_raw.rds")))
}

# qc table with explicit schema (format a)
qc_table <- tibble(sample=character(), before=integer(), after=integer(), removed=integer(), removed_pct=double())

raw_files <- list.files(out_load, "_raw.rds$", full.names = TRUE)

for (f in raw_files) {
  samp <- gsub("_raw.rds$", "", basename(f)); message("qc -> ", samp)
  so <- readRDS(f)

# remove empty barcodes
  if (!"nFeature_RNA" %in% colnames(so@meta.data)) so[["nFeature_RNA"]] <- Matrix::colSums(so@assays$RNA@counts > 0)
  so <- so[, so$nFeature_RNA > 0]

  if (is.null(so$percent.mt) || all(is.na(so$percent.mt))) so$percent.mt <- PercentageFeatureSet(so, pattern = "^mt-|^MT-|^mt")

  before_n <- ncol(so); message("  cells before qc: ", before_n)

  tryCatch({
    p1 <- VlnPlot(so, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3) + theme_bw() + ggtitle(paste0(samp, " before QC"))
    ggsave(file.path(out_qc_plots, paste0(samp,"_QC_before.png")), p1, width = 10, height = 4, dpi = 300)
  }, error = function(e) message("  warn: before-QC plot failed: ", e$message))

  so_f <- subset(so, subset =
                   nFeature_RNA >= params$min_features &
                   nCount_RNA >= params$min_counts &
                   percent.mt <= params$percent_mt_max)

  after_n <- ncol(so_f); message("  cells after qc: ", after_n)
  saveRDS(so_f, file.path(out_filtered, paste0(samp, "_filtered.rds")))
  message("  saved filtered rds: ", file.path(out_filtered, paste0(samp, "_filtered.rds")))

  tryCatch({
    p2 <- VlnPlot(so_f, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3) + theme_bw() + ggtitle(paste0(samp, " after QC"))
    ggsave(file.path(out_qc_plots, paste0(samp,"_QC_after.png")), p2, width = 10, height = 4, dpi = 300)
  }, error = function(e) message("  warn: after-QC plot failed: ", e$message))

  qc_table <- add_row(qc_table,
                      sample = samp,
                      before = before_n,
                      after = after_n,
                      removed = before_n - after_n,
                      removed_pct = round((before_n - after_n)/max(1,before_n)*100,2))
}

write_csv(as_tibble(params), file.path(out_qc, "parameters_stage1.csv"))
write_csv(qc_table, file.path(out_qc, "cell_counts_before_after_stage1.csv"))
message("pcos stage1 completed -> ", out_filtered)
