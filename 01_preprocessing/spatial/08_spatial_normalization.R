#!/usr/bin/env Rscript
# qc filtering + sctransform. pcos fix (ensembl->symbol merge) embedded.

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(Matrix); library(readr); library(tibble)
  library(viridisLite)
})

# config
data_root <- "~/projects/spatial/data"
pcos_data_dir <- file.path(data_root, "pcos_st")
aging_data_dir <- file.path(data_root, "aging_st")

load_root <- "~/projects/spatial/output/01_load"
filter_root <- "~/projects/spatial/output/02_qc"
sct_out_root <- "~/projects/spatial/output/03_sct"
pcos_fixed_dir <- file.path(sct_out_root, "pcos_fixed")
dir.create(filter_root, recursive=TRUE, showWarnings=FALSE)
dir.create(sct_out_root, recursive=TRUE, showWarnings=FALSE)
dir.create(pcos_fixed_dir, recursive=TRUE, showWarnings=FALSE)

# filtering thresholds (tweak if needed)
min_features <- 200
max_mt <- 20

# helper to get counts robustly
get_counts_from_so <- function(so) {
  cnt <- NULL
  try({ cnt <- GetAssayData(so, assay="RNA", layer="counts") }, silent=TRUE)
  if (is.null(cnt)) {
    if (!is.null(so@assays$RNA) && !is.null(so@assays$RNA@layers$counts)) cnt <- so@assays$RNA@layers$counts
  }
  if (is.null(cnt)) stop("Cannot find counts layer in Seurat object")
  if (!inherits(cnt, "dgCMatrix")) cnt <- Matrix(cnt, sparse=TRUE)
  return(cnt)
}

# pcos-specific function: map ens->symbol and merge duplicates
pcos_map_and_merge <- function(counts, feat_file) {
  feat <- read.table(feat_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
  colnames(feat)[1:2] <- c("ensembl","symbol")
  ens <- rownames(counts)
  map <- setNames(feat$symbol, feat$ensembl)[ens]
  map[is.na(map) | map==""] <- ens[is.na(map) | map==""]
  groups <- split(seq_along(map), map)
  merge_one <- function(idx) {
    if (length(idx) == 1) {
      m <- counts[idx, , drop=FALSE]
      if (!inherits(m, "dgCMatrix")) m <- Matrix(m, sparse=TRUE)
      return(m)
    } else {
      summed <- Matrix::colSums(counts[idx, , drop=FALSE])
      Matrix(summed, nrow=1, sparse=TRUE)
    }
  }
  merged_list <- lapply(groups, merge_one)
  merged_mat <- do.call(rbind, merged_list)
  rownames(merged_mat) <- names(groups)
  colnames(merged_mat) <- colnames(counts)
  return(merged_mat)
}

# process a sample: filters + sctransform (with option pcos_fix_flag)
process_sample <- function(raw_so, sample_name, sample_dir, pcos_fix=FALSE) {
# basic qc metrics already present
  raw_so[["percent.mt"]] <- PercentageFeatureSet(raw_so, pattern="^mt-|^MT-|^mt\\.")
  meta_before <- raw_so@meta.data %>% rownames_to_column("barcode")
# filtering cells
  keep_cells <- WhichCells(raw_so, expression = nFeature_RNA > min_features & percent.mt < max_mt)
  so_filt <- subset(raw_so, cells = keep_cells)
# save filtered object
  out_dir <- file.path(filter_root, sample_name, "filtered")
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  saveRDS(so_filt, file.path(out_dir, paste0(sample_name, "_filtered.rds")))
# qc summary csv
  qc_before <- tibble(sample=sample_name, total_cells=ncol(raw_so), filtered_cells=ncol(so_filt),
                      mean_nFeature_before=mean(raw_so$nFeature_RNA), mean_nFeature_after=mean(so_filt$nFeature_RNA))
  write_csv(qc_before, file.path(out_dir, paste0(sample_name, "_qc_summary_before_after.csv")))

# sctransform: for pcos we may need to rebuild counts with symbols
  if (pcos_fix) {
    feat_file <- file.path(sample_dir, "features.tsv.gz")
    counts <- get_counts_from_so(so_filt)
    merged <- pcos_map_and_merge(counts, feat_file)
# create new seurat object for sct
    so2 <- CreateSeuratObject(counts = merged, meta.data = so_filt@meta.data)
# preserve any spatial images
    if (length(so_filt@images) > 0) so2@images <- so_filt@images
    message("  PCOS fix: rebuilt counts to SYMBOLs, genes=", nrow(merged))
    so2 <- SCTransform(so2, vst.flavor="v2", verbose=FALSE, conserve.memory=TRUE)
    out_file <- file.path(pcos_fixed_dir, paste0(sample_name, "_filtered_sct_fixed.rds"))
    saveRDS(so2, out_file)
  } else {
# standard path
    so2 <- SCTransform(so_filt, vst.flavor="v2", verbose=FALSE, conserve.memory=TRUE)
    out_file <- file.path(sct_out_root, sample_name, paste0(sample_name, "_filtered_sct.rds"))
    dir.create(dirname(out_file), recursive=TRUE, showWarnings=FALSE)
    saveRDS(so2, out_file)
  }
  return(list(filtered = so_filt, sct = so2))
}

# run pcos samples
pcos_samples <- list.dirs(pcos_data_dir, recursive = FALSE, full.names = FALSE)
for (s in pcos_samples) {
  sd <- file.path(pcos_data_dir, s)
# load raw object from stage1 load dir
  raw_rds <- file.path("~/projects/spatial/output/01_load/pcos", s, paste0(s, "_raw.rds"))
  if (!file.exists(raw_rds)) {
    message("Missing stage1 raw rds:", raw_rds); next
  }
  raw_so <- readRDS(raw_rds)
  message("Processing PCOS sample: ", s)
  res <- process_sample(raw_so, s, sd, pcos_fix = TRUE)
}

# run aging samples (no special fix)
aging_samples <- list.dirs(aging_data_dir, recursive = FALSE, full.names = FALSE)
for (s in aging_samples) {
  sd <- file.path(aging_data_dir, s)
  raw_rds <- file.path("~/projects/spatial/output/01_load/aging", s, paste0(s, "_raw.rds"))
  if (!file.exists(raw_rds)) {
    message("Missing stage1 raw rds:", raw_rds); next
  }
  raw_so <- readRDS(raw_rds)
  message("Processing Aging sample: ", s)
  res <- process_sample(raw_so, s, sd, pcos_fix = FALSE)
}

message("Stage 2 complete: filtered and SCT objects saved.")