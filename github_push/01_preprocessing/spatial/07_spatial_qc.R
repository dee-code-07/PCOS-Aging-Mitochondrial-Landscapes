#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(png)
  library(tibble)
})

# user config
datasets <- list(
  pcos  = "~/projects/spatial/data/pcos_st",
  aging = "~/projects/spatial/data/aging_st"
)

out_base <- "~/projects/spatial/output"
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

min_features <- 100
min_counts   <- 500
max_mt       <- 20

qc_summary_csv <- file.path(out_base, "02_qc", "spatial_qc_summary_stage1_v4.csv")

# helpers

read_spatial_counts <- function(sdir) {
  h5f <- list.files(sdir, pattern = "filtered_feature_bc_matrix.*\\.h5$", full.names = TRUE)
  if (length(h5f) > 0) {
    message("    Using H5: ", basename(h5f[1]))
    raw <- Read10X_h5(h5f[1])
    if (is.list(raw)) {
      if ("Gene Expression" %in% names(raw)) return(raw[["Gene Expression"]])
      if ("RNA" %in% names(raw)) return(raw[["RNA"]])
      return(raw[[1]])
    } else return(raw)
  }
  mfile <- file.path(sdir, "matrix.mtx.gz")
  if (file.exists(mfile)) {
    message("    Using matrix.mtx.gz")
    m <- readMM(mfile)
    feats <- read.table(file.path(sdir, "features.tsv.gz"), sep = "\t", stringsAsFactors = FALSE)
    bcs   <- read.table(file.path(sdir, "barcodes.tsv.gz"), sep = "\t", stringsAsFactors = FALSE)
    rownames(m) <- feats$V1
    colnames(m) <- bcs$V1
    return(m)
  }
  stop("No valid count file found in: ", sdir)
}

read_feature_table <- function(sdir) {
  fts <- file.path(sdir, "features.tsv.gz")
  if (!file.exists(fts)) return(NULL)
  df <- read.table(fts, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  if (ncol(df) >= 2) {
    return(tibble(ensembl = df[,1], symbol = df[,2]))
  }
  return(NULL)
}

read_visium_positions <- function(spdir) {
  posf <- file.path(spdir, "tissue_positions_list.csv")
  if (!file.exists(posf)) return(NULL)
  df <- tryCatch(read.table(posf, sep = ",", header = FALSE, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || ncol(df) < 6) return(NULL)
  colnames(df)[1:6] <- c("barcode","in_tissue","array_row","array_col","pxl_col","pxl_row")
  return(df)
}

# robust percent.mt logic:
# - if rownames look like ensembl and symbol_map exists -> map and detect '^mt'
# - else if rownames have '^mt' directly -> use percentagefeatureset on rownames
# - else fallback to manual detection on symbol_map or rownames
compute_percent_mt_autodetect <- function(so, symbol_map = NULL) {
# get counts layer
  counts_mat <- tryCatch(GetAssayData(so, layer = "counts"), error = function(e) as.matrix(so@assays$RNA@counts))
  rnames <- rownames(counts_mat)

  looks_ensembl <- any(grepl("^ENSMUS", rnames, ignore.case = FALSE))
# case 1: rownames are symbols themselves and contain mt- => use percentagefeatureset on rownames
  if (any(grepl("^(mt|MT)", rnames))) {
# percentagefeatureset expects features in rownames
    pct <- tryCatch(PercentageFeatureSet(so, pattern = "^(mt|MT)"), error = function(e) {
# manual compute if necessary
      mt_idx <- grepl("^(mt|MT)", rnames)
      mt_counts <- Matrix::colSums(counts_mat[mt_idx, , drop = FALSE])
      tot_counts <- Matrix::colSums(counts_mat)
      100 * mt_counts / pmax(1, tot_counts)
    })
    return(pct)
  }

# case 2: rownames look like ensembl and symbol_map provided
  if (looks_ensembl && !is.null(symbol_map)) {
# map ensembl -> symbol
    syms <- symbol_map[rnames]
    syms[is.na(syms)] <- ""
    mt_idx <- grepl("^(mt|MT)", syms)
    if (any(mt_idx)) {
      mt_counts <- Matrix::colSums(counts_mat[mt_idx, , drop = FALSE])
      tot_counts <- Matrix::colSums(counts_mat)
      return(100 * mt_counts / pmax(1, tot_counts))
    }
  }

# case 3: symbol_map exists and contains mt entries in values
  if (!is.null(symbol_map)) {
# symbol_map is named vector ensembl->symbol
    syms_all <- symbol_map[rownames(counts_mat)]
    syms_all[is.na(syms_all)] <- ""
    if (any(grepl("^(mt|MT)", syms_all))) {
      mt_idx <- grepl("^(mt|MT)", syms_all)
      mt_counts <- Matrix::colSums(counts_mat[mt_idx, , drop = FALSE])
      tot_counts <- Matrix::colSums(counts_mat)
      return(100 * mt_counts / pmax(1, tot_counts))
    }
  }

# last resort: try matching common mitochondrial gene names (nd, co, cytb, atp) in symbols if available
  if (!is.null(symbol_map)) {
    syms_l <- tolower(symbol_map[rownames(counts_mat)])
    mito_pattern <- paste0("\\b(mt[-_]?nd|mt[-_]?co|mt[-_]?cytb|mt[-_]?atp|nd[1-6]|co[1-3]|cytb|atp6)\\b")
    mt_idx2 <- grepl(mito_pattern, syms_l)
    if (any(mt_idx2)) {
      mt_counts <- Matrix::colSums(counts_mat[mt_idx2, , drop = FALSE])
      tot_counts <- Matrix::colSums(counts_mat)
      return(100 * mt_counts / pmax(1, tot_counts))
    }
  }

# if nothing found, return zeros (safe)
  return(rep(0, ncol(counts_mat)))
}

# dot-only plot
plot_spots_only <- function(meta_df, expr_vec, title="", out_file=NULL) {
  df <- meta_df
  df$expr <- expr_vec[match(df$barcode, names(expr_vec))]
  df <- df[!is.na(df$x) & !is.na(df$y), ]
  g <- ggplot(df, aes(x=x, y=y, color=expr)) +
    geom_point(size=1.2) +
    scale_color_viridis_c(option = "C") +
    theme_void() +
    coord_equal() +
    ggtitle(title)
  if (!is.null(out_file)) ggsave(out_file, g, width = 6, height = 6, dpi = 300)
}

# main process
process_dataset <- function(name, base_dir, qc_list) {
  message("\n=== PROCESSING dataset: ", name, " ===\n")

  sample_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  if (length(sample_dirs) == 0) stop("No sample dirs in ", base_dir)

  out_load <- file.path(out_base, "01_load", name); dir.create(out_load, recursive = TRUE, showWarnings = FALSE)
  out_qc   <- file.path(out_base, "02_qc", name); dir.create(out_qc, recursive = TRUE, showWarnings = FALSE)
  qc_plots <- file.path(out_qc, "plots"); dir.create(qc_plots, recursive = TRUE, showWarnings = FALSE)
  qc_filt  <- file.path(out_qc, "filtered"); dir.create(qc_filt, recursive = TRUE, showWarnings = FALSE)
  sp_out   <- file.path(out_base, "05_spatial", name); dir.create(sp_out, recursive = TRUE, showWarnings = FALSE)

  for (sdir in sample_dirs) {
    samp <- basename(sdir)
    message(" - Sample: ", samp)

# load counts
    counts <- read_spatial_counts(sdir)

# create seurat object
    so <- CreateSeuratObject(counts = counts, project = name)
    so$orig.ident <- samp

# save raw loaded object (fix you asked)
    saveRDS(so, file.path(out_load, paste0(samp, "_raw.rds")))

# load feature table if exists (ensembl->symbol)
    ft <- read_feature_table(sdir)
    symbol_map <- NULL
    if (!is.null(ft)) symbol_map <- setNames(ft$symbol, ft$ensembl)

# compute qc metrics
    counts_mat <- GetAssayData(so, layer = "counts")
    so$nFeature_RNA <- Matrix::colSums(counts_mat > 0)
    so$nCount_RNA   <- Matrix::colSums(counts_mat)

# compute percent.mt robustly
    so$percent.mt <- compute_percent_mt_autodetect(so, symbol_map)

# before stats
    before_n <- ncol(so)

# attach spatial coords (if present)
    spdir <- file.path(sdir, "spatial")
    coords_df <- NULL
    has_visium_coords <- FALSE
    if (dir.exists(spdir)) {
      coords_df <- read_visium_positions(spdir)
      if (!is.null(coords_df) && nrow(coords_df) > 0) {
        coords_df <- coords_df[coords_df$in_tissue == 1, ]
        rownames(coords_df) <- coords_df$barcode
        so@meta.data$barcode <- colnames(so)
        so@meta.data$pxl_col <- coords_df[so@meta.data$barcode, "pxl_col"]
        so@meta.data$pxl_row <- coords_df[so@meta.data$barcode, "pxl_row"]
        has_visium_coords <- TRUE
      }
    }

# filtering
    keep_idx <- which(so$nFeature_RNA >= min_features & so$nCount_RNA >= min_counts & so$percent.mt <= max_mt)
    if (length(keep_idx) == 0) {
      message("   WARNING: no spots passed QC for ", samp, " — keeping all")
      so_f <- so
    } else {
      so_f <- subset(so, cells = colnames(so)[keep_idx])
    }

    after_n <- ncol(so_f)
    removed <- before_n - after_n
    removed_pct <- round(removed / max(1, before_n) * 100, 2)

# collect stats (some basic)
    stats <- tibble(
      dataset = name,
      sample = samp,
      before_n = before_n,
      after_n = after_n,
      removed = removed,
      removed_pct = removed_pct,
      before_median_nFeature = median(so$nFeature_RNA, na.rm = TRUE),
      after_median_nFeature = median(so_f$nFeature_RNA, na.rm = TRUE),
      before_median_nCount = median(so$nCount_RNA, na.rm = TRUE),
      after_median_nCount = median(so_f$nCount_RNA, na.rm = TRUE),
      before_median_percent_mt = median(so$percent.mt, na.rm = TRUE),
      after_median_percent_mt = median(so_f$percent.mt, na.rm = TRUE),
      before_mean_nFeature = mean(so$nFeature_RNA, na.rm = TRUE),
      after_mean_nFeature = mean(so_f$nFeature_RNA, na.rm = TRUE),
      before_mean_nCount = mean(so$nCount_RNA, na.rm = TRUE),
      after_mean_nCount = mean(so_f$nCount_RNA, na.rm = TRUE),
      before_mean_percent_mt = mean(so$percent.mt, na.rm = TRUE),
      after_mean_percent_mt = mean(so_f$percent.mt, na.rm = TRUE),
      before_min_percent_mt = min(so$percent.mt, na.rm = TRUE),
      before_max_percent_mt = max(so$percent.mt, na.rm = TRUE),
      after_min_percent_mt = min(so_f$percent.mt, na.rm = TRUE),
      after_max_percent_mt = max(so_f$percent.mt, na.rm = TRUE)
    )

    qc_list[[length(qc_list) + 1]] <- stats

# save filtered object
    saveRDS(so_f, file.path(qc_filt, paste0(samp, "_filtered.rds")))

# qc plots
    p_before <- VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
      ggtitle(paste0(samp, " - before QC")) + theme_bw()
    ggsave(file.path(qc_plots, paste0(samp, "_beforeQC.png")), p_before, width = 10, height = 4, dpi = 300)

    p_after <- VlnPlot(so_f, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
      ggtitle(paste0(samp, " - after QC")) + theme_bw()
    ggsave(file.path(qc_plots, paste0(samp, "_afterQC.png")), p_after, width = 10, height = 4, dpi = 300)

# dot-only spatial preview
    meta_df <- so_f@meta.data
    meta_df$barcode <- rownames(meta_df)
    if (has_visium_coords && !all(is.na(meta_df$pxl_col))) {
      meta_df$x <- as.numeric(meta_df$pxl_col)
      meta_df$y <- as.numeric(meta_df$pxl_row)
    } else {
      n <- nrow(meta_df)
      meta_df$x <- seq_len(n)
      meta_df$y <- rep(1, n)
    }
    expr_vec <- so_f$nCount_RNA
    names(expr_vec) <- colnames(so_f)
    out_png <- file.path(sp_out, paste0(samp, "_spots_only.png"))
    tryCatch({
      plot_spots_only(meta_df, expr_vec, title = paste0(samp, " - Spots Only"), out_file = out_png)
    }, error = function(e) message("   Spot-only plot failed for ", samp, ": ", e$message))

    message(" - Done sample: ", samp, " (kept: ", after_n, "/", before_n, ")\n")
  } # end samples

  return(qc_list)
}

# run
qc_list <- list()
for (nm in names(datasets)) qc_list <- process_dataset(nm, datasets[[nm]], qc_list)

qc_df <- bind_rows(qc_list)
dir.create(dirname(qc_summary_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(qc_df, qc_summary_csv, row.names = FALSE)
message("QC summary written to: ", qc_summary_csv)
message("STAGE-1 COMPLETE (v4).")
