#!/usr/bin/env Rscript
# task a v3: fixed spatial label transfer
# project: mitochondrial dysfunction in pcos and ovarian aging
# script:  a_stage5_label_transfer_fixed_v3.r
# fixes over v2:
# 1. hardcoded cell type column = "final_celltype" (confirmed from metadata)
# previously the auto-detector fell back to seurat_clusters (integers).
# 2. also checks aging annotated object for the same column name and warns
# if it differs.
# 3. duplicate sample guard: c2/c3 and p2/p3 share the same filtered rds —
# the script detects and reports this so you can decide whether to
# treat them as one sample or re-split upstream.
# outputs (per sample):
# - <sample>_with_celltypes.rds
# - <sample>_metadata_with_celltypes.csv
# - <sample>_celltype_summary.csv
# - <sample>_prediction_scores.csv
# - <sample>_spatial_celltype.png / .tiff
# - <sample>_prediction_score_dist.png
# summary:
# - label_transfer_summary.csv
# - prediction_score_overview.png / .tiff

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(patchwork)
})

# configuration

PROJECT_ROOT <- "E:/Documents/mini_project"

ref_pcos_path  <- file.path(PROJECT_ROOT, "scrna/output/07_annotation/pcos/pcos_annotated.rds")
ref_aging_path <- file.path(PROJECT_ROOT, "scrna/output/07_annotation/aging/aging_annotated.rds")

spatial_qc_base   <- file.path(PROJECT_ROOT, "spatial/output/02_qc")
spatial_raw_pcos  <- file.path(PROJECT_ROOT, "spatial/raw/GSE296728_PCOS_ST")

out_base <- file.path(PROJECT_ROOT, "spatial/output/06_label_transfer_FIXED")
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

SCORE_THRESHOLD <- 0.5

# cell type column (hardcoded after inspection)
# confirmed from: colnames(ref@meta.data) — the annotated column is
# "final_celltype" in both pcos and aging objects.
# change these if your aging object uses a different column name.
CELLTYPE_COL_PCOS  <- "final_celltype"
CELLTYPE_COL_AGING <- "final_celltype"   # will be verified at runtime

# explicit mapping: sample id → raw folder → tar.gz file prefix (inside the tar)
# the tar.gz contains: barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
pcos_sample_map <- list(
  c1 = list(folder = "pcost_ctrl_1",   tar = "GSM8975360_C1_matrix.tar.gz",   sym = "C1"),
  c2 = list(folder = "pcost_ctrl_2_3", tar = "GSM8975361_C2_matrix.tar.gz",   sym = "C2"),
  c3 = list(folder = "pcost_ctrl_2_3", tar = "GSM8975361_C3_matrix.tar.gz",   sym = "C3"),
  p1 = list(folder = "pcost_case_1",   tar = "GSM8975362_P1_matrix.tar.gz",   sym = "P1"),
  p2 = list(folder = "pcost_case_2_3", tar = "GSM8975363_P2_matrix.tar.gz",   sym = "P2"),
  p3 = list(folder = "pcost_case_2_3", tar = "GSM8975363_P3_matrix.tar.gz",   sym = "P3")
)

samples <- list(
  pcos  = c("c1", "c2", "c3", "p1", "p2", "p3"),
  aging = c("ya_1", "ya_2", "ya_3", "ya_4")
)

# helper: verify celltype column exists

verify_celltype_col <- function(so, col_name, dataset_name) {
  if (!col_name %in% colnames(so@meta.data)) {
    message("\n  ERROR: Column '", col_name, "' not found in ", dataset_name)
    message("  Available columns: ", paste(colnames(so@meta.data), collapse = ", "))
    stop("Cell type column '", col_name, "' missing from ", dataset_name,
         ". Update CELLTYPE_COL_PCOS / CELLTYPE_COL_AGING in CONFIG.")
  }
  vals <- sort(unique(so@meta.data[[col_name]]))
  message("  [OK] '", col_name, "' found in ", dataset_name,
          " (", length(vals), " cell types)")
  message("  Cell types: ", paste(vals, collapse = ", "))
}

# helper: extract features.tsv.gz from tar.gz
# the pcos spatial raw data stores the matrix in .tar.gz archives.
# we extract features.tsv.gz to a temp dir and read it from there.

extract_features_from_tar <- function(tar_path, temp_dir) {
  
  if (!file.exists(tar_path)) {
    message("  WARNING: tar.gz not found: ", tar_path)
    return(NULL)
  }
  
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  
# list contents of the tar to find the features file
  listing <- tryCatch(
    system2("tar", args = c("-tzf", shQuote(tar_path)), stdout = TRUE, stderr = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(listing) || length(listing) == 0) {
# windows fallback: use r's untar
    tryCatch(
      untar(tar_path, exdir = temp_dir),
      error = function(e) message("  WARNING: untar failed: ", e$message)
    )
  } else {
# extract only the features file
    feat_entry <- listing[grepl("features", listing, ignore.case = TRUE)]
    if (length(feat_entry) == 0) feat_entry <- listing  # extract all if not found
    tryCatch(
      untar(tar_path, files = feat_entry, exdir = temp_dir),
      error = function(e) {
        message("  Selective extract failed, trying full extract...")
        untar(tar_path, exdir = temp_dir)
      }
    )
  }
  
# find the extracted features file
  feat_file <- list.files(temp_dir, pattern = "features", recursive = TRUE,
                          full.names = TRUE)
  if (length(feat_file) == 0) {
    message("  WARNING: features file not found after extraction in ", temp_dir)
    return(NULL)
  }
  
  message("  Extracted: ", basename(feat_file[1]))
  return(feat_file[1])
}

# helper: map ensembl → symbol

aggregate_rows_sparse <- function(counts, newrownames) {
  if (length(unique(newrownames)) == length(newrownames)) {
    rownames(counts) <- newrownames
    return(counts)
  }
  sm    <- summary(counts)
  uniq  <- unique(newrownames)
  new_i <- match(newrownames[sm$i], uniq)
  sparseMatrix(
    i = new_i, j = sm$j, x = sm$x,
    dims = c(length(uniq), ncol(counts)),
    dimnames = list(uniq, colnames(counts))
  )
}

map_pcos_genes <- function(sp, s) {
  
  info     <- pcos_sample_map[[s]]
  tar_path <- file.path(spatial_raw_pcos, info$folder, info$tar)
  temp_dir <- file.path(tempdir(), paste0("pcos_feat_", s))
  
  message("  Looking for tar: ", basename(tar_path))
  feat_file <- extract_features_from_tar(tar_path, temp_dir)
  
  if (is.null(feat_file)) {
    message("  Cannot map genes for ", s, " — skipping (rownames kept as-is)")
    return(sp)
  }
  
# read features table
  feat <- tryCatch(
    read.delim(feat_file, header = FALSE, stringsAsFactors = FALSE),
    error = function(e) {
      read.csv(feat_file, header = FALSE, stringsAsFactors = FALSE)
    }
  )
  
  if (ncol(feat) < 2) {
    message("  WARNING: features file has unexpected format")
    return(sp)
  }
  colnames(feat)[1:2] <- c("ensembl", "symbol")
  
# get counts from spatial object
  DefaultAssay(sp) <- "RNA"
  counts <- tryCatch(
    GetAssayData(sp, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(sp, assay = "RNA", slot = "counts")
  )
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  
  old    <- rownames(counts)
  mapped <- feat$symbol[match(old, feat$ensembl)]
# fallback: keep original if no mapping found
  mapped[is.na(mapped) | mapped == ""] <- old[is.na(mapped) | mapped == ""]
  
  n_mapped <- sum(mapped != old)
  message("  Mapped ", n_mapped, "/", length(old), " ENSEMBL IDs to gene symbols")
  
  agg <- aggregate_rows_sparse(counts, mapped)
  sp[["RNA"]] <- CreateAssayObject(counts = agg)
  
# clean up temp
  unlink(temp_dir, recursive = TRUE)
  
  return(sp)
}

# helper: check gene overlap before anchoring

check_gene_overlap <- function(ref, sp, min_overlap = 200) {
  ref_genes <- rownames(ref)
  sp_genes  <- rownames(sp)
  overlap   <- intersect(ref_genes, sp_genes)
  message("  Gene overlap: ", length(overlap), " shared genes",
          " (ref=", length(ref_genes), ", query=", length(sp_genes), ")")
  if (length(overlap) < min_overlap) {
    message("  WARNING: Only ", length(overlap), " shared genes — transfer may fail.")
    message("  First 10 query genes: ", paste(head(sp_genes, 10), collapse = ", "))
    message("  First 10 ref genes:   ", paste(head(ref_genes, 10), collapse = ", "))
    return(FALSE)
  }
  return(TRUE)
}

# helper: safe metadata to csv
# avoids the "barcode column already exists" error

safe_meta_to_csv <- function(so, filepath) {
  meta <- so@meta.data
# only add barcode column if it doesn't already exist
  if (!"barcode" %in% colnames(meta)) {
    meta <- rownames_to_column(meta, "barcode")
  }
  write_csv(meta, filepath)
}

# helper: prepare reference

prep_reference <- function(ref, celltype_col) {
  DefaultAssay(ref) <- "RNA"
# collapse seurat v5 layers if needed
  if (inherits(ref[["RNA"]], "Assay5")) {
    tryCatch(ref <- JoinLayers(ref, assay = "RNA"), error = function(e) NULL)
  }
  ref <- NormalizeData(ref, verbose = FALSE)
  ref <- FindVariableFeatures(ref, nfeatures = 3000, verbose = FALSE)
  ref <- ScaleData(ref, verbose = FALSE)
  ref <- RunPCA(ref, verbose = FALSE)
  ref$celltype_for_transfer <- ref@meta.data[[celltype_col]]
  Idents(ref) <- ref$celltype_for_transfer
  return(ref)
}

prep_spatial_query <- function(sp) {
  DefaultAssay(sp) <- "RNA"
  sp <- NormalizeData(sp, verbose = FALSE)
  sp <- FindVariableFeatures(sp, nfeatures = 3000, verbose = FALSE)
  sp <- ScaleData(sp, verbose = FALSE)
  sp <- RunPCA(sp, verbose = FALSE)
  return(sp)
}

# load and prepare references

message("\n", strrep("=", 70))
message("  LOADING REFERENCE OBJECTS")
message(strrep("=", 70))

message("\nLoading PCOS reference...")
ref_pcos <- readRDS(ref_pcos_path)
verify_celltype_col(ref_pcos, CELLTYPE_COL_PCOS, "PCOS")

message("\nLoading Aging reference...")
ref_aging <- readRDS(ref_aging_path)
# verify aging — fall back gracefully if column name differs
if (!CELLTYPE_COL_AGING %in% colnames(ref_aging@meta.data)) {
  message("  WARNING: '", CELLTYPE_COL_AGING, "' not found in Aging object.")
  message("  Available: ", paste(colnames(ref_aging@meta.data), collapse = ", "))
# try to find it automatically
  char_cols <- colnames(ref_aging@meta.data)[
    sapply(colnames(ref_aging@meta.data), function(x)
      is.character(ref_aging@meta.data[[x]]) || is.factor(ref_aging@meta.data[[x]]))]
  n_uniq <- sapply(char_cols, function(x) length(unique(ref_aging@meta.data[[x]])))
  reasonable <- char_cols[n_uniq >= 3 & n_uniq <= 50]
  if (length(reasonable) > 0) {
    CELLTYPE_COL_AGING <- reasonable[which.max(n_uniq[reasonable])]
    message("  Auto-selected aging column: '", CELLTYPE_COL_AGING, "'")
  } else {
    stop("Cannot find cell type column in Aging object. Check colnames(ref_aging@meta.data).")
  }
}
verify_celltype_col(ref_aging, CELLTYPE_COL_AGING, "Aging")

message("\nPreparing references...")
ref_pcos  <- prep_reference(ref_pcos,  CELLTYPE_COL_PCOS)
ref_aging <- prep_reference(ref_aging, CELLTYPE_COL_AGING)

# main loop

summary_rows <- list()

for (dataset in names(samples)) {
  
  message("\n", strrep("█", 70))
  message("  DATASET: ", toupper(dataset))
  message(strrep("█", 70))
  
  ref        <- if (dataset == "pcos") ref_pcos  else ref_aging
  ct_col_ref <- if (dataset == "pcos") CELLTYPE_COL_PCOS else CELLTYPE_COL_AGING
  
  for (s in samples[[dataset]]) {
    
    message("\n--- Sample: ", s, " ---")
    
# load filtered spatial object
    filt_dir <- file.path(spatial_qc_base, dataset, "filtered")
    fpath    <- list.files(filt_dir, pattern = paste0("^", s, "_filtered\\.rds$"),
                           full.names = TRUE)
    if (length(fpath) == 0) {
# try looser match
      fpath <- list.files(filt_dir, pattern = s, full.names = TRUE)
      fpath <- fpath[grep("\\.rds$", fpath)]
    }
    if (length(fpath) == 0) {
      message("  SKIP: RDS not found in ", filt_dir); next
    }
    
    sp <- readRDS(fpath[1])
    message("  Loaded: ", basename(fpath[1]), " — ", ncol(sp), " spots")
    
# duplicate sample guard
# c2/c3 share pcost_ctrl_2_3; p2/p3 share pcost_case_2_3.
# if this sample's rds is identical to a previously processed one,
# warn the user — results will be duplicated.
    rds_key <- normalizePath(fpath[1], mustWork = FALSE)
    if (exists("processed_rds_paths") && rds_key %in% processed_rds_paths) {
      message("  ⚠ WARNING: '", basename(fpath[1]), "' was already used for a previous",
              " sample. Results for '", s, "' will be identical to that sample.")
      message("  This means c2/c3 or p2/p3 were never split into separate RDS files.")
      message("  Proceeding — but consider re-checking your spatial QC stage.")
    }
    if (!exists("processed_rds_paths")) processed_rds_paths <- character(0)
    processed_rds_paths <- c(processed_rds_paths, rds_key)
    
# add x/y coordinates for plotting
    if ("pxl_col" %in% colnames(sp@meta.data)) {
      sp@meta.data$x <- as.numeric(sp@meta.data$pxl_col)
      sp@meta.data$y <- as.numeric(sp@meta.data$pxl_row)
    }
    
# pcos: map ensembl → gene symbols
    if (dataset == "pcos") {
      message("  Mapping ENSEMBL → gene symbols...")
      sp <- map_pcos_genes(sp, s)
    }
    
# check gene overlap before running pca
    ok <- check_gene_overlap(ref, sp)
    if (!ok) {
      message("  SKIP: insufficient gene overlap for ", s)
      next
    }
    
# normalize and reduce spatial query
    message("  Normalizing spatial query...")
    sp <- prep_spatial_query(sp)
    
# findtransferanchors
    message("  Finding transfer anchors...")
    anchors <- tryCatch(
      FindTransferAnchors(
        reference            = ref,
        query                = sp,
        normalization.method = "LogNormalize",
        dims                 = 1:30,
        verbose              = FALSE
      ),
      error = function(e) { message("  ERROR: ", e$message); NULL }
    )
    if (is.null(anchors)) { message("  SKIP: anchor finding failed"); next }
    
# transferdata
    message("  Transferring labels...")
    preds <- tryCatch(
      TransferData(
        anchorset = anchors,
        refdata   = ref$celltype_for_transfer,
        dims      = 1:30,
        verbose   = FALSE
      ),
      error = function(e) { message("  ERROR: ", e$message); NULL }
    )
    if (is.null(preds)) { message("  SKIP: transfer failed"); next }
    
    sp$predicted_celltype <- preds$predicted.id
    sp$prediction_score   <- preds$prediction.score.max
    
# qc metrics
    n_total  <- ncol(sp)
    n_pass   <- sum(sp$prediction_score >= SCORE_THRESHOLD, na.rm = TRUE)
    pct_pass <- round(100 * n_pass / n_total, 1)
    med_sc   <- round(median(sp$prediction_score, na.rm = TRUE), 3)
    message(sprintf("  QC: %d/%d spots pass (%.1f%%) | median score = %.3f",
                    n_pass, n_total, pct_pass, med_sc))
    
# save outputs
    out_dir <- file.path(out_base, dataset, s)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
# rds
    saveRDS(sp, file.path(out_dir, paste0(s, "_with_celltypes.rds")))
    
# metadata csv — safe version avoids barcode duplication
    safe_meta_to_csv(sp, file.path(out_dir, paste0(s, "_metadata_with_celltypes.csv")))
    
# cell type frequency table
    ct_tab <- sp@meta.data %>%
      count(predicted_celltype) %>%
      mutate(fraction = round(n / sum(n), 4)) %>%
      arrange(desc(n))
    write_csv(ct_tab, file.path(out_dir, paste0(s, "_celltype_summary.csv")))
    
# per-cell-type prediction score summary
    score_tab <- sp@meta.data %>%
      group_by(predicted_celltype) %>%
      summarise(
        n_spots       = n(),
        median_score  = round(median(prediction_score, na.rm = TRUE), 4),
        mean_score    = round(mean(prediction_score,   na.rm = TRUE), 4),
        pct_above_0.5 = round(100 * mean(prediction_score >= 0.5, na.rm = TRUE), 1),
        .groups = "drop"
      ) %>% arrange(desc(median_score))
    write_csv(score_tab, file.path(out_dir, paste0(s, "_prediction_scores.csv")))
    
# plot 1: spatial cell type map
    if (all(c("x","y") %in% colnames(sp@meta.data))) {
      
      n_ct    <- length(unique(sp$predicted_celltype))
      ct_pal  <- if (n_ct <= 8)  RColorBrewer::brewer.pal(max(3, n_ct), "Set2") else
        if (n_ct <= 12) RColorBrewer::brewer.pal(n_ct, "Set3") else
          scales::hue_pal()(n_ct)
      
      p_spatial <- ggplot(sp@meta.data,
                          aes(x = x, y = -y, colour = predicted_celltype)) +
        geom_point(size = 0.7, alpha = 0.85) +
        scale_colour_manual(values = ct_pal, name = "Cell Type") +
        coord_fixed() +
        theme_void(base_size = 10) +
        theme(
          legend.position  = "right",
          legend.text      = element_text(size = 7),
          legend.title     = element_text(size = 8, face = "bold"),
          plot.title       = element_text(size = 11, face = "bold", hjust = 0.5),
          plot.subtitle    = element_text(size = 8,  hjust = 0.5, colour = "grey50"),
          plot.margin      = margin(5, 5, 5, 5)
        ) +
        labs(
          title    = paste0(toupper(dataset), " | ", s, "  —  Predicted Cell Types"),
          subtitle = sprintf("%d/%d spots score ≥ %.2f  |  median = %.3f",
                             n_pass, n_total, SCORE_THRESHOLD, med_sc)
        )
      
      ggsave(file.path(out_dir, paste0(s, "_spatial_celltype.png")),
             p_spatial, width = 8, height = 7, dpi = 300)
      ggsave(file.path(out_dir, paste0(s, "_spatial_celltype.tiff")),
             p_spatial, width = 8, height = 7, dpi = 600, compression = "lzw")
    }
    
# plot 2: prediction score violin by cell type
    p_scores <- tryCatch({
      ggplot(sp@meta.data,
             aes(x = reorder(predicted_celltype, prediction_score, FUN = median),
                 y = prediction_score, fill = predicted_celltype)) +
        geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +
        geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.7) +
        geom_hline(yintercept = SCORE_THRESHOLD, linetype = "dashed",
                   colour = "firebrick", linewidth = 0.6) +
        coord_flip() +
        theme_classic(base_size = 10) +
        theme(legend.position = "none",
              axis.text.y     = element_text(size = 8),
              plot.title      = element_text(face = "bold", hjust = 0.5)) +
        labs(title   = paste0(s, "  —  Prediction Score by Cell Type"),
             x = NULL, y = "Max Prediction Score",
             caption = paste0("Red dashed line = threshold (", SCORE_THRESHOLD, ")"))
    }, error = function(e) NULL)
    
    if (!is.null(p_scores))
      ggsave(file.path(out_dir, paste0(s, "_prediction_score_dist.png")),
             p_scores, width = 7, height = 5, dpi = 300)
    
# summary row
    summary_rows[[length(summary_rows) + 1]] <- tibble(
      dataset          = dataset,
      sample           = s,
      n_spots_total    = n_total,
      n_spots_pass     = n_pass,
      pct_pass         = pct_pass,
      median_score     = med_sc,
      n_celltypes      = length(unique(sp$predicted_celltype)),
      top_celltype     = ct_tab$predicted_celltype[1],
      top_celltype_pct = round(100 * ct_tab$n[1] / n_total, 1)
    )
    
    message("  ✓ Done: ", out_dir)
    gc()
  }
}

# summary

if (length(summary_rows) == 0) {
  message("\nNo samples completed successfully. Check errors above.")
  quit(status = 1)
}

summary_df <- bind_rows(summary_rows)
write_csv(summary_df, file.path(out_base, "label_transfer_summary.csv"))

message("\n", strrep("=", 70))
message("  LABEL TRANSFER SUMMARY")
message(strrep("=", 70))
print(as.data.frame(summary_df))

# cross-sample prediction score overview
all_scores <- list()
for (dataset in names(samples)) {
  for (s in samples[[dataset]]) {
    rds_path <- file.path(out_base, dataset, s, paste0(s, "_with_celltypes.rds"))
    if (!file.exists(rds_path)) next
    sp_tmp <- readRDS(rds_path)
    all_scores[[paste0(dataset, "_", s)]] <- data.frame(
      sample   = paste0(dataset, "_", s),
      dataset  = dataset,
      score    = sp_tmp$prediction_score,
      celltype = sp_tmp$predicted_celltype,
      stringsAsFactors = FALSE
    )
    rm(sp_tmp); gc()
  }
}

if (length(all_scores) > 0) {
  all_scores_df <- bind_rows(all_scores)
  
  p_overview <- ggplot(all_scores_df,
                       aes(x = sample, y = score, fill = dataset)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.4) +
    geom_hline(yintercept = SCORE_THRESHOLD, linetype = "dashed",
               colour = "firebrick", linewidth = 0.7) +
    scale_fill_manual(values = c(pcos = "#D55E00", aging = "#0072B2")) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      plot.title      = element_text(face = "bold", hjust = 0.5),
      legend.position = "top"
    ) +
    labs(
      title   = "Prediction Score Distribution — All Samples",
      x = NULL, y = "Max Prediction Score",
      caption = paste0("Red dashed = threshold (", SCORE_THRESHOLD, ")")
    )
  
  ggsave(file.path(out_base, "prediction_score_overview.png"),
         p_overview, width = 10, height = 5, dpi = 300)
  ggsave(file.path(out_base, "prediction_score_overview.tiff"),
         p_overview, width = 10, height = 5, dpi = 600, compression = "lzw")
  
  message("✓ Overview plot saved.")
}

message("\n✓ Summary: ", file.path(out_base, "label_transfer_summary.csv"))
message("✓ TASK A COMPLETE. Outputs: ", out_base, "\n")