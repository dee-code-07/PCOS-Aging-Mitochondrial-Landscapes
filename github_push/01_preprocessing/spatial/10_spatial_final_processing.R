#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(Matrix)
})

# user options
base_dir <- path.expand("~/projects/spatial")
integration_dir <- file.path(base_dir, "output", "04_integration")
filtered_base <- file.path(base_dir, "output", "02_qc")
sct_base <- file.path(base_dir, "output", "03_sct")
data_base <- file.path(base_dir, "data")
out_base <- file.path(base_dir, "output", "05_spatial_visualization_final_with_clusters")
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# markers (symbol list) - edit if needed
markers <- c(
  "TOMM20","HSPD1","NDUFA9","COX4I1",
  "SOD2","PRDX3","PINK1","BNIP3","SQSTM1",
  "CDKN1A","CDKN2A","TP53","SERPINE1",
  "IL6","TNF","IL1B","CCL2","CXCL1","CXCL2"
)

datasets <- list(
  pcos = c("c1","c2","c3","p1","p2","p3"),
  aging = c("ya_1","ya_2","ya_3","ya_4")
)

# helpers
read_features_map <- function(fpath) {
  if (!file.exists(fpath)) return(NULL)
  df <- tryCatch(read.table(fpath, sep = "\t", header = FALSE, stringsAsFactors = FALSE, quote="", comment.char=""), error = function(e) NULL)
  if (is.null(df) || ncol(df) < 2) return(NULL)
  setNames(df[,2], df[,1])  # names = ensembl, values = symbol
}

get_matrix_for_sample <- function(sct_obj, spatial_obj) {
  m <- NULL
  if (!is.null(sct_obj)) {
    if ("SCT" %in% Assays(sct_obj)) m <- tryCatch(GetAssayData(sct_obj, assay="SCT", slot="data"), error=function(e) NULL)
    if (is.null(m) && "RNA" %in% Assays(sct_obj)) m <- tryCatch(GetAssayData(sct_obj, assay="RNA", slot="data"), error=function(e) NULL)
  }
  if (is.null(m) && !is.null(spatial_obj) && "RNA" %in% Assays(spatial_obj)) {
    m <- tryCatch(GetAssayData(spatial_obj, assay="RNA", slot="counts"), error=function(e) NULL)
    if (!is.null(m)) m <- Matrix::log1p(m)
  }
  return(m)
}

plot_spots_expr <- function(meta, expr_vec, title_txt, outfile, width=6, height=6, dpi=600) {
  d <- meta
  d$expr <- expr_vec[match(d$barcode, names(expr_vec))]
  d$expr[is.na(d$expr)] <- 0
  p <- ggplot(d, aes(x=x, y=y, color = expr)) +
    geom_point(size = 0.9) +
    scale_color_viridis_c(option="C") +
    theme_void() + coord_equal() +
    ggtitle(title_txt) + theme(plot.title = element_text(hjust = 0.5), legend.position="right")
  ggsave(outfile, p, width = width, height = height, dpi = dpi)
}

# robust barcode matcher: returns named vector cluster_for_barcodes (names = spatial barcodes)
map_clusters_to_barcodes <- function(int_obj, spatial_barcodes) {
# int_obj colnames are cell ids (e.g. "aaaca...-1_1"); int_obj$seurat_clusters exists
  int_cells <- colnames(int_obj)
  clusters <- as.character(int_obj$seurat_clusters)
  names(clusters) <- int_cells

  out <- setNames(rep(NA_character_, length(spatial_barcodes)), spatial_barcodes)

# prepare helper variants
  for (i in seq_along(spatial_barcodes)) {
    b <- spatial_barcodes[i]
    candidates <- unique(c(
      b,
      sub("-1$","", b),
      paste0(b, "_1"),
      paste0(sub("-1$","", b), "_1"),
      paste0(b, "_2"), # sometimes sample index differs; try _2 as fallback
      paste0(sub("-1$","", b), "_2"),
      sub("_1$", "", int_cells) # include int variants fallback
    ))
    found <- NA
    for (cand in candidates) {
      if (is.na(cand) || cand == "") next
      if (cand %in% int_cells) { found <- clusters[cand]; break }
    }
# last resort: match by partial barcode (without suffix) to first matching int_cell
    if (is.na(found)) {
      base_strip <- sub("-1$","", b)
      idx <- which(grepl(base_strip, int_cells, fixed=TRUE))
      if (length(idx) > 0) found <- clusters[int_cells[idx[1]]]
    }
    out[i] <- ifelse(is.null(found) || length(found)==0, NA_character_, as.character(found))
  }
  return(out)
}

# main
for (ds in names(datasets)) {
  int_rds <- file.path(integration_dir, ds, paste0(ds, "_integrated.rds"))
  if (!file.exists(int_rds)) {
    message("Integrated object not found for ", ds, " -> skipping.")
    next
  }
  message("Loading integrated object for: ", ds)
  int_obj <- readRDS(int_rds)

  for (samp in datasets[[ds]]) {
    message("\n=== Sample: ", samp, " (", ds, ") ===")
# spatial filtered rds
    filt_dir <- file.path(filtered_base, ds, "filtered")
    filt_files <- list.files(filt_dir, pattern = samp, full.names = TRUE)
    if (length(filt_files) == 0) {
      message("  Filtered spatial not found for ", samp, " -> skipping.")
      next
    }
    spatial_rds <- filt_files[1]
    spatial_obj <- readRDS(spatial_rds)

# load sct if exists
    sct_dir <- file.path(sct_base, ds)
    sct_files <- list.files(sct_dir, pattern = samp, full.names = TRUE)
    sct_obj <- if (length(sct_files)>0) readRDS(sct_files[1]) else NULL

# ensure meta coords present
    meta <- spatial_obj@meta.data
    if ("pxl_col" %in% colnames(meta) && "pxl_row" %in% colnames(meta)) {
      meta$x <- as.numeric(meta$pxl_col); meta$y <- as.numeric(meta$pxl_row)
    } else if ("x" %in% colnames(meta) && "y" %in% colnames(meta)) {
      meta$x <- as.numeric(meta$x); meta$y <- as.numeric(meta$y)
    } else {
      meta$x <- seq_len(nrow(meta)); meta$y <- rep(1, nrow(meta))
    }
    if (!("barcode" %in% colnames(meta))) meta$barcode <- rownames(meta)
    spatial_obj@meta.data <- meta

# map clusters: very robust mapping function
    mapped_clusters <- map_clusters_to_barcodes(int_obj, colnames(spatial_obj))
    spatial_obj@meta.data$cluster_mapped <- mapped_clusters[match(spatial_obj@meta.data$barcode, names(mapped_clusters))]
    spatial_obj@meta.data$cluster_mapped[is.na(spatial_obj@meta.data$cluster_mapped)] <- "unassigned"

# write updated annotated rds
    out_sample_dir <- file.path(out_base, ds, samp)
    dir.create(out_sample_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(spatial_obj, file.path(out_sample_dir, paste0(samp, "_spatial_with_clusters.rds")))

# save metadata csv
    write.csv(spatial_obj@meta.data, file.path(out_sample_dir, paste0(samp, "_spots_metadata_with_clusters.csv")), row.names = FALSE)

# read feature map for sample (pcos has features.tsv.gz)
    feat_file <- file.path(data_base, ifelse(ds=="pcos", "pcos_st", "aging_st"), samp, "features.tsv.gz")
    feat_map <- if (file.exists(feat_file)) read.table(feat_file, sep="\t", header=FALSE, stringsAsFactors=FALSE, quote="", comment.char="") else NULL
    feat_named <- NULL
    if (!is.null(feat_map) && ncol(feat_map) >= 2) feat_named <- setNames(feat_map[,2], feat_map[,1])

# expression matrix to use
    mat <- get_matrix_for_sample(sct_obj, spatial_obj)
    if (is.null(mat)) {
      message("  No expression matrix found for ", samp, " -> skipping marker plots.")
      next
    }

# per-gene plotting + summaries
    combined_df <- NULL
    for (g in markers) {
      expr_vec <- NULL
# direct symbol in mat
      if (g %in% rownames(mat)) {
        v <- as.numeric(mat[g,]); names(v) <- colnames(mat); expr_vec <- v
      } else if (!is.null(feat_named)) {
        ens_hits <- names(feat_named)[tolower(feat_named) == tolower(g)]
        if (length(ens_hits) > 0) {
          rows_pres <- intersect(ens_hits, rownames(mat))
          if (length(rows_pres) > 0) {
            v <- Matrix::colSums(as.matrix(mat[rows_pres, , drop=FALSE])); names(v) <- colnames(mat); expr_vec <- v
          }
        }
      }
# strip versions
      if (is.null(expr_vec)) {
        rn_strip <- sub("\\.\\d+$","", rownames(mat))
        idx <- which(tolower(rn_strip) == tolower(g))
        if (length(idx)>0) { v <- Matrix::colSums(as.matrix(mat[idx, , drop=FALSE])); names(v) <- colnames(mat); expr_vec <- v }
      }
      if (is.null(expr_vec)) {
        message("  Marker not found for sample ", samp, ": ", g); next
      }

# map expr to spatial barcodes
      expr_spatial <- rep(NA_real_, ncol(spatial_obj))
      names(expr_spatial) <- colnames(spatial_obj)
      common <- intersect(names(expr_vec), names(expr_spatial))
      if (length(common)>0) expr_spatial[common] <- expr_vec[common]
      else {
        ss_n <- sub("-1$","", names(expr_vec)); ss_sp <- sub("-1$","", names(expr_spatial))
        for (i in seq_along(ss_sp)) {
          idx <- which(ss_n == ss_sp[i])
          if (length(idx)>0) expr_spatial[i] <- expr_vec[idx[1]]
        }
      }

# plot
      meta_df <- spatial_obj@meta.data
      pngf <- file.path(out_sample_dir, paste0(samp, "_", g, "_spots.png"))
      try(plot_spots_expr(meta_df, expr_spatial, paste0(samp, " : ", g), pngf), silent=TRUE)

# cluster summary
      sumdf <- meta_df %>%
        mutate(expr = expr_spatial[match(barcode, names(expr_spatial))]) %>%
        group_by(cluster_mapped) %>%
        summarise(mean_expr = mean(expr, na.rm = TRUE),
                  detected_frac = sum(!is.na(expr) & expr > 0)/n(),
                  n_spots = n()) %>%
        arrange(desc(mean_expr)) %>%
        as.data.frame()
      write.csv(sumdf, file.path(out_sample_dir, paste0(samp, "_", g, "_summary_by_cluster.csv")), row.names = FALSE)

# accumulate combined
      df2 <- sumdf[, c("cluster_mapped", "mean_expr")]
      colnames(df2) <- c("cluster", g)
      if (is.null(combined_df)) combined_df <- df2 else combined_df <- merge(combined_df, df2, by="cluster", all=TRUE)
    }

    if (!is.null(combined_df)) write.csv(combined_df, file.path(out_sample_dir, paste0(samp, "_combined_marker_means_by_cluster.csv")), row.names = FALSE)

    message("  Finished sample: ", samp, " -> outputs in: ", out_sample_dir)
  } # end samples
} # end datasets

message("All done. Outputs: ", out_base)