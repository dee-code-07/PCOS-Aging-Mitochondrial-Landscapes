#!/usr/bin/env Rscript
# task c v3: spatially variable genes — manual moran's i
# project: mitochondrial dysfunction in pcos and ovarian aging
# script:  c_spatially_variable_genes_v3.r
# why v1 and v2 failed:
# your spatial pipeline uses createseuratobject() + coordinate metadata,
# not load10x_spatial(). this means no @images slot ever exists.
# findspatiallyvariablefeatures() requires a proper visium image object
# and silently returns 0 genes without one.
# solution:
# implement moran's i from scratch using:
# 1. pxl_col / pxl_row coordinates from object metadata
# 2. a k-nearest-neighbour spatial weights matrix (k=6, visium neighbours)
# 3. gene expression from the sct assay
# this is exactly what seurat does internally — we just do it directly.
# performance:
# using sparse matrix operations. expect ~2-5 min per sample for 2000 genes.
# inputs:
# spatial/output/03_sct/pcos/*.rds       (per-sample sct objects)
# spatial/output/03_sct/aging/*.rds
# analysis/09_gene_prioritization/top_priority_genes_tier1.csv
# outputs (analysis/16_spatially_variable_genes/):
# pcos_svgs_per_sample.csv
# pcos_svgs_top50.csv
# aging_svgs_per_sample.csv
# aging_svgs_top50.csv
# svg_priority_overlap.csv
# pcos_svg_spatialfeatureplot.png/.tiff
# aging_svg_spatialfeatureplot.png/.tiff
# svg_overlap_barplot.png/.tiff
# svg_moran_dotplot.png/.tiff
# svg_run_summary.csv

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(patchwork)
  library(RColorBrewer)
})

# configuration

PROJECT_ROOT  <- "E:/Documents/mini_project"

SCT_PCOS_DIR  <- file.path(PROJECT_ROOT, "spatial/output/03_sct/pcos")
SCT_AGING_DIR <- file.path(PROJECT_ROOT, "spatial/output/03_sct/aging")

PRIORITY_PATH <- file.path(PROJECT_ROOT,
                           "analysis/09_gene_prioritization/Top_Priority_Genes_Tier1.csv")

OUT_DIR <- file.path(PROJECT_ROOT, "analysis/16_spatially_variable_genes")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_GENES_TEST <- 2000   # top variable genes to test per sample
N_SVG_REPORT <- 50     # genes in output tables
N_SVG_PLOT   <- 9      # genes in spatial plots (3×3)
KNN_K        <- 6      # spatial neighbours (6 = Visium hex grid)
N_CORES      <- 1      # keep at 1 for Windows (no fork())

set.seed(42)

# theme

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = base_size + 1),
      plot.subtitle = element_text(hjust = 0.5, colour = "grey40", size = base_size - 1),
      axis.text     = element_text(colour = "black"),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text    = element_text(face = "bold")
    )
}

# core: build knn spatial weight matrix
# given spot coordinates, build a row-normalised k-nn weight matrix w.
# w[i,j] = 1/k if j is among k nearest neighbours of i, else 0.
# this is the standard spatial weights matrix for moran's i.

build_spatial_weights <- function(x, y, k = 6) {
  
  n   <- length(x)
  coords <- cbind(x, y)
  
  message("  Building ", k, "-NN spatial weights for ", n, " spots...")
  
# compute pairwise distances efficiently using vectorised ops
# for large n, use a chunked approach to avoid n×n matrix in ram
# here we build a sparse knn weight matrix directly
  
# use kd-tree via built-in dist if n is manageable, else chunk
  if (n <= 5000) {
    D <- as.matrix(dist(coords, method = "euclidean"))
# for each row, find k nearest (excluding self)
    W_i <- integer(0)
    W_j <- integer(0)
    for (i in seq_len(n)) {
      nbrs <- order(D[i, ])[2:(k + 1)]  # exclude self (rank 1)
      W_i  <- c(W_i, rep(i, k))
      W_j  <- c(W_j, nbrs)
    }
  } else {
# chunked for larger datasets
    message("  Large dataset — using chunked kNN (", n, " spots)...")
    chunk_size <- 500
    W_i <- integer(0)
    W_j <- integer(0)
    for (start in seq(1, n, by = chunk_size)) {
      end   <- min(start + chunk_size - 1, n)
      chunk <- coords[start:end, , drop = FALSE]
      D_chunk <- as.matrix(dist(rbind(chunk, coords)))[
        seq_len(end - start + 1), , drop = FALSE]
      for (ii in seq_len(nrow(D_chunk))) {
        i    <- start + ii - 1
        nbrs <- order(D_chunk[ii, ])[2:(k + 1)]
        W_i  <- c(W_i, rep(i, k))
        W_j  <- c(W_j, nbrs)
      }
    }
  }
  
# build sparse weight matrix, row-normalised (each row sums to 1)
  W <- sparseMatrix(
    i    = W_i,
    j    = W_j,
    x    = rep(1 / k, length(W_i)),
    dims = c(n, n)
  )
  
  message("  Weight matrix built: ", n, "×", n, " sparse")
  return(W)
}

# core: compute moran's i for multiple genes
# vectorised moran's i:
# i = (n / s0) * (z'wz / z'z)
# where z = gene expression - mean, s0 = sum of all weights
# returns a data.frame with gene, morans_i, expected_i, z_score, p_value

compute_morans_i_batch <- function(expr_mat, W, gene_names = NULL) {
# expr_mat: genes × spots (dense or sparse)
# w:        spots × spots sparse weight matrix (row-normalised)
  
  if (is.null(gene_names)) gene_names <- rownames(expr_mat)
  
  n  <- ncol(expr_mat)   # number of spots
  S0 <- sum(W)           # sum of all weights
  
# convert to dense for vectorised z-score computation
# process in chunks to avoid memory blow-up
  chunk_size <- 200
  n_genes    <- nrow(expr_mat)
  
  results <- vector("list", ceiling(n_genes / chunk_size))
  chunk_idx <- 0
  
  for (start in seq(1, n_genes, by = chunk_size)) {
    chunk_idx <- chunk_idx + 1
    end       <- min(start + chunk_size - 1, n_genes)
    
# extract chunk as dense matrix: genes × spots
    E <- as.matrix(expr_mat[start:end, , drop = FALSE])
    
# centre each gene (subtract mean)
    gene_means <- rowMeans(E)
    Z <- E - gene_means   # genes × spots
    
# compute z'z (variance term) per gene
    zz <- rowSums(Z^2)
    
# compute z'wz (spatial autocorrelation term) per gene
# z %*% t(w) gives genes × spots; then row-wise dot with z
# = rowsums(z * (z %*% t(w)))   but w is spots×spots
# z'wz for one gene = z %*% w %*% z
# vectorised: (z %*% w) * z, sum rows
    ZW   <- Z %*% W           # genes × spots  (dense × sparse = dense)
    zWz  <- rowSums(ZW * Z)   # genes
    
# moran's i
    I_obs <- (n / S0) * (zWz / pmax(zz, .Machine$double.eps))
    
# expected value under null: e[i] = -1/(n-1)
    E_I <- -1 / (n - 1)
    
# variance of i under normality assumption (cliff & ord 1981)
# simplified formula: var(i) ≈ (n^2*s1 - n*s2 + 3*s0^2) / (s0^2*(n^2-1))
# where s1 = 0.5 * sum((w + w')^2), s2 = sum((rowsums(w)+colsums(w))^2)
# for row-normalised w with k=6: s0=n, so use fast approximation
    S1    <- 0.5 * sum((W + t(W))^2)
    rW    <- Matrix::rowSums(W) + Matrix::colSums(W)
    S2    <- sum(rW^2)
    Var_I <- (n^2 * S1 - n * S2 + 3 * S0^2) /
      (S0^2 * (n^2 - 1)) -
      E_I^2
    
    SD_I  <- sqrt(max(Var_I, .Machine$double.eps))
    z_scores <- (I_obs - E_I) / SD_I
    p_vals   <- 2 * pnorm(-abs(z_scores))   # two-tailed
    
    results[[chunk_idx]] <- data.frame(
      gene      = gene_names[start:end],
      morans_i  = round(I_obs, 6),
      expected_i = round(E_I, 6),
      z_score   = round(z_scores, 4),
      p_value   = p_vals,
      stringsAsFactors = FALSE
    )
    
    if (chunk_idx %% 5 == 0)
      message("    ... ", end, "/", n_genes, " genes processed")
  }
  
  result_df <- bind_rows(results)
  return(result_df)
}

# helper: find and load sct files

find_sct_files <- function(sct_dir) {
  patterns <- c("_filtered_sct_fixed\\.rds$", "_filtered_sct\\.rds$",
                "_sct_fixed\\.rds$", "_sct\\.rds$")
  for (pat in patterns) {
    files <- list.files(sct_dir, pattern = pat, full.names = TRUE)
    files <- files[!grepl("backup|BACKUP|old", basename(files))]
    if (length(files) > 0) return(files)
  }
# last resort
  list.files(sct_dir, pattern = "\\.rds$", full.names = TRUE)
}

# helper: extract coordinates from object

get_spatial_coords <- function(so, sample_name) {
  
  meta <- so@meta.data
  
# try coordinate column pairs in priority order
  coord_pairs <- list(
    c("pxl_col", "pxl_row"),
    c("x", "y"),
    c("array_col", "array_row"),
    c("imagerow", "imagecol"),
    c("row", "col")
  )
  
  for (pair in coord_pairs) {
    cx <- pair[1]; cy <- pair[2]
    if (all(c(cx, cy) %in% colnames(meta))) {
      x <- as.numeric(meta[[cx]])
      y <- as.numeric(meta[[cy]])
      if (!all(is.na(x)) && length(unique(x)) > 1) {
        message("  Coordinates from: '", cx, "' / '", cy, "'")
        return(list(x = x, y = y, xcol = cx, ycol = cy))
      }
    }
  }
  
# try @images slot as last resort
  if (length(Images(so)) > 0) {
    img_name <- Images(so)[1]
    img_obj  <- so@images[[img_name]]
    if (!is.null(img_obj@coordinates)) {
      coords <- img_obj@coordinates
      if (all(c("imagerow","imagecol") %in% colnames(coords))) {
        message("  Coordinates from @images slot")
        return(list(
          x    = coords$imagecol[match(rownames(meta), rownames(coords))],
          y    = coords$imagerow[match(rownames(meta), rownames(coords))],
          xcol = "imagecol", ycol = "imagerow"
        ))
      }
    }
  }
  
  message("  WARNING: No spatial coordinates found for ", sample_name)
  message("  Available metadata cols: ", paste(colnames(meta), collapse = ", "))
  return(NULL)
}

# main: run svg per sample

run_svg_dataset <- function(sct_dir, dataset_label) {
  
  message("\n", strrep("█", 65))
  message("  DATASET: ", dataset_label)
  message(strrep("█", 65))
  
  rds_files <- find_sct_files(sct_dir)
  message("  Found ", length(rds_files), " sample files:")
  for (f in rds_files) message("    ", basename(f))
  
  all_svg_results <- list()
  
  for (f in rds_files) {
    
# clean up sample name
    samp <- basename(f)
    samp <- sub("_filtered_sct_fixed\\.rds$|_filtered_sct\\.rds$|_sct_fixed\\.rds$|_sct\\.rds$", "", samp)
    samp <- sub("\\.rds$", "", samp)
    
    message("\n  ── Sample: ", samp, " ──")
    
    so <- tryCatch(readRDS(f),
                   error = function(e) { message("  Cannot load: ", e$message); NULL })
    if (is.null(so)) next
    
    message("  Loaded: ", ncol(so), " spots × ", nrow(so), " genes")
    
# get spatial coordinates
    coords <- get_spatial_coords(so, samp)
    if (is.null(coords)) {
      message("  SKIP — no coordinates found")
      next
    }
    
    x <- coords$x
    y <- coords$y
    
# remove spots with na coordinates
    valid <- !is.na(x) & !is.na(y)
    if (sum(valid) < 50) {
      message("  SKIP — fewer than 50 spots with valid coordinates (",
              sum(valid), ")")
      next
    }
    if (sum(!valid) > 0) {
      message("  Removing ", sum(!valid), " spots with NA coordinates")
      so <- so[, valid]
      x  <- x[valid]
      y  <- y[valid]
    }
    
    message("  Valid spots for SVG: ", length(x))
    
# determine assay and get expression
    use_assay <- NULL
    for (a in c("SCT", "Spatial", "RNA")) {
      if (a %in% names(so@assays)) { use_assay <- a; break }
    }
    DefaultAssay(so) <- use_assay
    message("  Assay: ", use_assay)
    
# collapse seurat v5 layers
    if (inherits(so[[use_assay]], "Assay5")) {
      tryCatch(so <- JoinLayers(so, assay = use_assay), error = function(e) NULL)
    }
    
# select genes to test
# use variable features if available, else top expressed genes
    vf <- VariableFeatures(so)
    if (length(vf) == 0) {
      so  <- FindVariableFeatures(so, nfeatures = N_GENES_TEST, verbose = FALSE)
      vf  <- VariableFeatures(so)
    }
    test_genes <- head(vf, N_GENES_TEST)
    message("  Testing ", length(test_genes), " variable genes")
    
# get normalised expression matrix (genes × spots)
    expr_mat <- tryCatch(
      GetAssayData(so, assay = use_assay, layer = "data")[test_genes, ],
      error = function(e)
        GetAssayData(so, assay = use_assay, slot  = "data")[test_genes, ]
    )
    
# convert to dense — needed for matrix ops
# only keep genes expressed in >5% of spots to reduce noise
    pct_expressed <- Matrix::rowMeans(expr_mat > 0)
    expr_mat <- expr_mat[pct_expressed > 0.05, ]
    message("  Genes after expression filter (>5% spots): ", nrow(expr_mat))
    
    if (nrow(expr_mat) == 0) {
      message("  SKIP — no genes passed expression filter")
      next
    }
    
# convert to dense matrix for moran's i computation
    expr_dense <- as.matrix(expr_mat)
    
# build spatial weights
    W <- tryCatch(
      build_spatial_weights(x, y, k = KNN_K),
      error = function(e) {
        message("  Weight matrix failed: ", e$message); NULL
      }
    )
    if (is.null(W)) next
    
# compute moran's i
    message("  Computing Moran's I for ", nrow(expr_dense), " genes...")
    t_start <- proc.time()
    
    moran_df <- tryCatch(
      compute_morans_i_batch(expr_dense, W),
      error = function(e) {
        message("  Moran's I computation failed: ", e$message); NULL
      }
    )
    
    t_elapsed <- round((proc.time() - t_start)["elapsed"], 1)
    message("  Completed in ", t_elapsed, "s")
    
    if (is.null(moran_df)) next
    
# multiple testing correction and ranking
    moran_df <- moran_df %>%
      mutate(
        p_adj = p.adjust(p_value, method = "BH"),
        spatially_variable = p_adj < 0.05 & morans_i > 0
      ) %>%
      arrange(desc(morans_i))
    
    n_sig <- sum(moran_df$spatially_variable, na.rm = TRUE)
    message("  Significant SVGs (FDR<0.05, I>0): ", n_sig)
    message("  Top 5 by Moran's I: ",
            paste(head(moran_df$gene[moran_df$morans_i > 0], 5), collapse = ", "))
    
# add sample info
    moran_df$sample  <- samp
    moran_df$dataset <- dataset_label
    moran_df$rank    <- seq_len(nrow(moran_df))
    
    all_svg_results[[samp]] <- list(
      moran_df = moran_df,
      coords   = data.frame(barcode = colnames(so), x = x, y = y),
      so       = so,
      assay    = use_assay
    )
    
    gc()
  }
  
  return(all_svg_results)
}

# aggregate across samples

aggregate_svg_results <- function(all_results, dataset_label) {
  
  if (length(all_results) == 0) return(NULL)
  
# only use samples that produced results
  valid <- Filter(function(r) !is.null(r$moran_df) && nrow(r$moran_df) > 0,
                  all_results)
  if (length(valid) == 0) return(NULL)
  
  message("\n  Aggregating across ", length(valid), " samples (", dataset_label, ")...")
  
# combine all per-sample tables
  per_sample_df <- bind_rows(lapply(valid, function(r) r$moran_df))
  
# aggregate: mean moran's i and count of samples where fdr < 0.05
  agg_df <- per_sample_df %>%
    group_by(gene) %>%
    summarise(
      n_samples_tested = n(),
      n_samples_sig    = sum(spatially_variable, na.rm = TRUE),
      mean_morans_i    = round(mean(morans_i, na.rm = TRUE), 5),
      max_morans_i     = round(max(morans_i,  na.rm = TRUE), 5),
      mean_z_score     = round(mean(z_score,  na.rm = TRUE), 3),
      min_p_adj        = round(min(p_adj,     na.rm = TRUE), 6),
      .groups = "drop"
    ) %>%
# rank: prioritise genes sig in most samples, then highest moran's i
    arrange(desc(n_samples_sig), desc(mean_morans_i)) %>%
    mutate(aggregated_rank = row_number())
  
  message("  Unique genes tested: ", nrow(agg_df))
  message("  SVGs significant in ≥1 sample: ",
          sum(agg_df$n_samples_sig >= 1))
  message("  SVGs significant in ALL samples: ",
          sum(agg_df$n_samples_sig == length(valid)))
  message("  Top 10 SVGs: ",
          paste(head(agg_df$gene, 10), collapse = ", "))
  
  return(list(
    per_sample = per_sample_df,
    aggregated = agg_df,
    top_genes  = head(agg_df$gene[agg_df$n_samples_sig >= 1], N_SVG_REPORT),
    all_results = valid
  ))
}

# run

message("\n", strrep("=", 65))
message("  TASK C v3: MANUAL MORAN'S I SVG ANALYSIS")
message(strrep("=", 65))

results_pcos  <- run_svg_dataset(SCT_PCOS_DIR,  "PCOS")
results_aging <- run_svg_dataset(SCT_AGING_DIR, "Aging")

agg_pcos  <- aggregate_svg_results(results_pcos,  "PCOS")
agg_aging <- aggregate_svg_results(results_aging, "Aging")

# save tables

message("\n", strrep("=", 65))
message("  SAVING RESULTS")
message(strrep("=", 65))

if (!is.null(agg_pcos)) {
  write_csv(agg_pcos$per_sample,
            file.path(OUT_DIR, "PCOS_SVGs_per_sample.csv"))
  write_csv(head(agg_pcos$aggregated, N_SVG_REPORT),
            file.path(OUT_DIR, "PCOS_SVGs_top50.csv"))
  message("✓ PCOS: ", nrow(agg_pcos$aggregated), " unique SVGs saved")
}

if (!is.null(agg_aging)) {
  write_csv(agg_aging$per_sample,
            file.path(OUT_DIR, "Aging_SVGs_per_sample.csv"))
  write_csv(head(agg_aging$aggregated, N_SVG_REPORT),
            file.path(OUT_DIR, "Aging_SVGs_top50.csv"))
  message("✓ Aging: ", nrow(agg_aging$aggregated), " unique SVGs saved")
}

# priority gene overlap

priority_genes <- tryCatch({
  df       <- read_csv(PRIORITY_PATH, show_col_types = FALSE)
  gene_col <- intersect(c("gene","Gene","gene_symbol","Gene_Symbol","symbol"),
                        colnames(df))[1]
  unique(df[[gene_col]])
}, error = function(e) { message("Priority genes: ", e$message); character(0) })

message("\n  Tier 1 priority genes loaded: ", length(priority_genes))

if (length(priority_genes) > 0) {
  
  pcos_genes  <- if (!is.null(agg_pcos))  agg_pcos$aggregated$gene  else character(0)
  aging_genes <- if (!is.null(agg_aging)) agg_aging$aggregated$gene else character(0)
  
# use all svgs for overlap (not just top 50)
  ov_pcos  <- intersect(priority_genes, pcos_genes)
  ov_aging <- intersect(priority_genes, aging_genes)
  ov_both  <- intersect(ov_pcos, ov_aging)
  
  message("  Tier1 ∩ PCOS SVGs:  ", length(ov_pcos))
  message("  Tier1 ∩ Aging SVGs: ", length(ov_aging))
  message("  Tier1 ∩ Both:       ", length(ov_both))
  if (length(ov_both) > 0)
    message("  Shared genes: ", paste(ov_both, collapse = ", "))
  
  all_hits <- union(ov_pcos, ov_aging)
  
  if (length(all_hits) > 0) {
    lookup_p <- agg_pcos$aggregated
    lookup_a <- agg_aging$aggregated
    
    overlap_df <- tibble(gene = all_hits) %>%
      left_join(lookup_p %>% select(gene, pcos_rank = aggregated_rank,
                                    pcos_morans_i = mean_morans_i,
                                    pcos_n_sig = n_samples_sig),
                by = "gene") %>%
      left_join(lookup_a %>% select(gene, aging_rank = aggregated_rank,
                                    aging_morans_i = mean_morans_i,
                                    aging_n_sig = n_samples_sig),
                by = "gene") %>%
      mutate(
        svg_in_pcos  = !is.na(pcos_rank),
        svg_in_aging = !is.na(aging_rank),
        svg_in_both  = gene %in% ov_both
      ) %>%
      arrange(pcos_rank)
    
    write_csv(overlap_df, file.path(OUT_DIR, "SVG_priority_overlap.csv"))
    message("✓ Overlap table: ", nrow(overlap_df), " genes")
  } else {
# still write an empty file so downstream scripts don't error
    write_csv(
      tibble(gene=character(), svg_in_pcos=logical(), svg_in_aging=logical(),
             svg_in_both=logical()),
      file.path(OUT_DIR, "SVG_priority_overlap.csv")
    )
    message("  No Tier1 genes found in SVG lists — empty overlap file written")
  }
  
# overlap barplot
  overlap_counts <- data.frame(
    category = c("SVG in PCOS only", "SVG in Aging only", "SVG in Both"),
    n        = c(length(setdiff(ov_pcos, ov_aging)),
                 length(setdiff(ov_aging, ov_pcos)),
                 length(ov_both))
  )
  
  p_overlap <- ggplot(overlap_counts,
                      aes(x = reorder(category, n), y = n, fill = category)) +
    geom_bar(stat = "identity", width = 0.55) +
    geom_text(aes(label = n), hjust = -0.3, fontface = "bold", size = 4.5) +
    scale_fill_manual(values = c("SVG in PCOS only"  = "#D55E00",
                                 "SVG in Aging only"  = "#0072B2",
                                 "SVG in Both"        = "#009E73")) +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    theme_pub() +
    theme(legend.position = "none") +
    labs(title    = "Tier 1 Priority Genes with Spatial Variability",
         subtitle = paste0("Total Tier 1 genes = ", length(priority_genes),
                           " | Method: Manual Moran's I"),
         x = NULL, y = "Number of genes")
  
  ggsave(file.path(OUT_DIR, "SVG_overlap_barplot.png"),
         p_overlap, width = 7, height = 4, dpi = 300)
  ggsave(file.path(OUT_DIR, "SVG_overlap_barplot.tiff"),
         p_overlap, width = 7, height = 4, dpi = 600, compression = "lzw")
  message("✓ Overlap barplot saved")
}

# spatial feature plots

plot_top_svgs_scatter <- function(agg_result, dataset_label) {
  
  if (is.null(agg_result) || length(agg_result$top_genes) == 0) return(invisible(NULL))
  
# find best sample (most svgs, largest)
  best_samp <- names(which.max(sapply(agg_result$all_results,
                                      function(r) sum(r$moran_df$spatially_variable))))
  if (is.null(best_samp)) best_samp <- names(agg_result$all_results)[1]
  
  r_best    <- agg_result$all_results[[best_samp]]
  so_best   <- r_best$so
  coords_df <- r_best$coords
  use_assay <- r_best$assay
  
  DefaultAssay(so_best) <- use_assay
  if (inherits(so_best[[use_assay]], "Assay5"))
    tryCatch(so_best <- JoinLayers(so_best, assay = use_assay), error = function(e) NULL)
  
# genes to plot: top significant svgs present in this sample
  plot_genes <- intersect(agg_result$top_genes, rownames(so_best))
  plot_genes <- head(plot_genes, N_SVG_PLOT)
  
  if (length(plot_genes) == 0) {
    message("  No plottable genes for ", dataset_label); return(invisible(NULL))
  }
  
  message("  Plotting ", length(plot_genes), " SVGs for ", dataset_label,
          " (sample: ", best_samp, ")")
  
# fetch expression
  expr_df <- tryCatch(
    as.data.frame(t(as.matrix(
      GetAssayData(so_best, assay = use_assay, layer = "data")[plot_genes, , drop = FALSE]
    ))),
    error = function(e)
      as.data.frame(t(as.matrix(
        GetAssayData(so_best, assay = use_assay, slot = "data")[plot_genes, , drop = FALSE]
      )))
  )
  expr_df$barcode <- rownames(expr_df)
  
# merge with coordinates
  plot_df <- left_join(coords_df, expr_df, by = "barcode")
  
  ncols  <- 3
  nrows  <- ceiling(length(plot_genes) / ncols)
  
  plot_list <- lapply(plot_genes, function(g) {
    ggplot(plot_df, aes(x = x, y = -y, colour = .data[[g]])) +
      geom_point(size = 0.7, alpha = 0.9) +
      scale_colour_gradientn(
        colours = c("lightgrey", "#fee090", "#fc8d59", "#d73027"),
        name    = "Expr", na.value = "lightgrey"
      ) +
      coord_fixed() +
      theme_void(base_size = 9) +
      theme(
        legend.key.height = unit(0.35, "cm"),
        legend.text       = element_text(size = 6.5),
        plot.title        = element_text(face = "bold.italic", size = 9, hjust = 0.5)
      ) +
      labs(title = g)
  })
  
  p_combined <- wrap_plots(plot_list, ncol = ncols) +
    plot_annotation(
      title    = paste0(dataset_label, "  —  Top Spatially Variable Genes"),
      subtitle = paste0("Sample: ", best_samp,
                        "  |  Method: Manual Moran's I (k=", KNN_K, " NN)"),
      theme    = theme(
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, colour = "grey50", size = 9)
      )
    )
  
  fname <- paste0(dataset_label, "_SVG_SpatialFeaturePlot")
  ggsave(file.path(OUT_DIR, paste0(fname, ".png")),
         p_combined, width = 4 * ncols, height = 4 * nrows + 0.5, dpi = 300)
  ggsave(file.path(OUT_DIR, paste0(fname, ".tiff")),
         p_combined, width = 4 * ncols, height = 4 * nrows + 0.5, dpi = 600,
         compression = "lzw")
  message("✓ Spatial feature plot saved: ", fname)
}

plot_top_svgs_scatter(agg_pcos,  "PCOS")
plot_top_svgs_scatter(agg_aging, "Aging")

# moran's i dotplot

make_moran_dotplot <- function(agg_pcos, agg_aging, n_show = 25) {
  
  build_df <- function(agg, label) {
    if (is.null(agg)) return(NULL)
    head(agg$aggregated, n_show) %>%
      mutate(dataset = label,
             gene    = factor(gene, levels = rev(gene)))
  }
  
  df_p <- build_df(agg_pcos,  "PCOS")
  df_a <- build_df(agg_aging, "Aging")
  
  make_panel <- function(df, colour, label) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    ggplot(df, aes(x = mean_morans_i, y = gene)) +
      geom_segment(aes(x = 0, xend = mean_morans_i, yend = gene),
                   colour = "grey82", linewidth = 0.5) +
      geom_point(aes(size = n_samples_sig, colour = mean_morans_i)) +
      scale_colour_gradient(low = "grey75", high = colour, guide = "none") +
      scale_size_continuous(name = "# Samples\n(FDR<0.05)",
                            range = c(2, 6), breaks = c(1, 2, 3, 4)) +
      geom_vline(xintercept = 0, linetype = "dashed",
                 colour = "grey60", linewidth = 0.4) +
      theme_pub(base_size = 9) +
      theme(axis.text.y = element_text(face = "italic", size = 7.5)) +
      labs(title = paste0(label, " — Top ", n_show, " SVGs"),
           x = "Mean Moran's I", y = NULL)
  }
  
  p1 <- make_panel(df_p, "#D55E00", "PCOS")
  p2 <- make_panel(df_a, "#0072B2", "Aging")
  plots <- Filter(Negate(is.null), list(p1, p2))
  if (length(plots) == 0) return(invisible(NULL))
  
  p_combined <- wrap_plots(plots, ncol = 2)
  ggsave(file.path(OUT_DIR, "SVG_moran_dotplot.png"),
         p_combined, width = 13, height = 8, dpi = 300)
  ggsave(file.path(OUT_DIR, "SVG_moran_dotplot.tiff"),
         p_combined, width = 13, height = 8, dpi = 600, compression = "lzw")
  message("✓ Moran's I dotplot saved")
}

make_moran_dotplot(agg_pcos, agg_aging)

# final summary

message("\n", strrep("=", 65))
message("  TASK C COMPLETE")
message(strrep("=", 65))

summary_df <- tibble(
  dataset           = c("PCOS", "Aging"),
  n_samples_run     = c(length(results_pcos), length(results_aging)),
  n_svgs_total      = c(if (!is.null(agg_pcos))  nrow(agg_pcos$aggregated)  else 0,
                        if (!is.null(agg_aging)) nrow(agg_aging$aggregated) else 0),
  n_svgs_sig_any    = c(if (!is.null(agg_pcos))
    sum(agg_pcos$aggregated$n_samples_sig >= 1)  else 0,
    if (!is.null(agg_aging))
      sum(agg_aging$aggregated$n_samples_sig >= 1) else 0),
  top3_svgs         = c(
    if (!is.null(agg_pcos)  && length(agg_pcos$top_genes)  > 0)
      paste(head(agg_pcos$top_genes,  3), collapse = ", ") else "none",
    if (!is.null(agg_aging) && length(agg_aging$top_genes) > 0)
      paste(head(agg_aging$top_genes, 3), collapse = ", ") else "none"
  )
)

print(as.data.frame(summary_df))
write_csv(summary_df, file.path(OUT_DIR, "SVG_run_summary.csv"))
message("\nOutputs in: ", OUT_DIR)