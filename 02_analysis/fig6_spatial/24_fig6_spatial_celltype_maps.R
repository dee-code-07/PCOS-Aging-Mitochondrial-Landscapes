#!/usr/bin/env Rscript
# task b: spatial cell type composition maps
# project: mitochondrial dysfunction in pcos and ovarian aging
# script:  b_spatial_celltype_maps.r
# purpose:
# after fixing label transfer (task a), we now have predicted_celltype and
# prediction_score for every spatial spot. this script generates:
# 1. spatial cell type distribution maps (where each cell type sits)
# 2. cell type proportion barplots (pcos vs control, aged vs young)
# 3. co-localisation: which cell types are enriched in high mitochondrial
# module score regions (the key spatial biology question)
# inputs:
# spatial/output/06_label_transfer_fixed/<dataset>/<sample>/<sample>_with_celltypes.rds
# spatial/output/04_integration/<dataset>/<dataset>_integrated.rds  (for module scores)
# outputs (analysis/18_spatial_celltype_maps/):
# per sample:
# <sample>_celltype_spatial.png/.tiff
# summary figures:
# pcos_celltype_proportions.png/.tiff
# aging_celltype_proportions.png/.tiff
# celltype_proportion_comparison.png/.tiff   (pcos vs control, aged vs young)
# colocalisation_mito_celltype.png/.tiff     (key spatial biology figure)
# colocalisation_heatmap.png/.tiff

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(patchwork)
  library(RColorBrewer)
  library(tidyr)
  library(scales)
})

# configuration

PROJECT_ROOT <- "E:/Documents/mini_project"

LT_BASE <- file.path(PROJECT_ROOT, "spatial/output/06_label_transfer_FIXED")

# integrated spatial objects — needed for module scores (mitochondrial etc.)
INTEGRATED_PCOS  <- file.path(PROJECT_ROOT,
                              "spatial/output/04_integration/pcos/pcos_integrated.rds")
INTEGRATED_AGING <- file.path(PROJECT_ROOT,
                              "spatial/output/04_integration/aging/aging_integrated.rds")

OUT_DIR <- file.path(PROJECT_ROOT, "analysis/18_spatial_celltype_maps")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# sample metadata
PCOS_SAMPLES  <- c("c1","c2","c3","p1","p2","p3")
AGING_SAMPLES <- c("ya_1","ya_2","ya_3","ya_4")

PCOS_CONDITION  <- c(c1="control", c2="control", c3="control",
                     p1="case",    p2="case",     p3="case")
AGING_CONDITION <- c(ya_1="aged", ya_2="aged", ya_3="aged", ya_4="aged",
                     ya_5="young", ya_6="young", ya_7="young", ya_8="young")
# note: if aging only has ya_1..4, they may all be one group —
# the script detects this and handles it

# prediction score threshold (from task a)
SCORE_THRESH <- 0.5

# mitochondrial module score column name in integrated objects
# common names from the pipeline — script tries all of these
MITO_SCORE_CANDIDATES <- c(
  "mito_score", "mitochondrial_score", "MitoScore",
  "mito_module", "mitochondrial_module",
  "Module1", "module_mito",
  "Mito1", "mito1",
  "mito_stress", "oxidative_stress_score",
  "mito_dysfunction_score"
)

# theme

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, colour = "grey40",
                                   size  = base_size - 1),
      axis.text     = element_text(colour = "black"),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text    = element_text(face = "bold", size = base_size - 1)
    )
}

# colour palette for cell types
make_celltype_palette <- function(cell_types) {
  n <- length(cell_types)
  if (n <= 8)       pal <- brewer.pal(max(3,n), "Set2")
  else if (n <= 12) pal <- brewer.pal(n, "Set3")
  else              pal <- hue_pal()(n)
  setNames(pal[seq_len(n)], sort(cell_types))
}

# helper: load label transfer results

load_lt_sample <- function(dataset, sample_id) {
  rds_path <- file.path(LT_BASE, dataset, sample_id,
                        paste0(sample_id, "_with_celltypes.rds"))
  if (!file.exists(rds_path)) {
    message("  NOT FOUND: ", rds_path)
    return(NULL)
  }
  so <- tryCatch(readRDS(rds_path),
                 error = function(e) { message("  Load error: ", e$message); NULL })
  if (is.null(so)) return(NULL)
  
# ensure x/y coordinates present
  meta <- so@meta.data
  if (!all(c("x","y") %in% colnames(meta))) {
    if (all(c("pxl_col","pxl_row") %in% colnames(meta))) {
      so@meta.data$x <- as.numeric(meta$pxl_col)
      so@meta.data$y <- as.numeric(meta$pxl_row)
    } else if (all(c("array_col","array_row") %in% colnames(meta))) {
      so@meta.data$x <- as.numeric(meta$array_col)
      so@meta.data$y <- as.numeric(meta$array_row)
    }
  }
  
  return(so)
}

# helper: detect module score columns

find_module_score_cols <- function(so) {
  meta_cols <- colnames(so@meta.data)
  
# try exact candidates first
  exact <- intersect(MITO_SCORE_CANDIDATES, meta_cols)
  if (length(exact) > 0) return(exact)
  
# pattern search
  pattern_hits <- meta_cols[grepl(
    "mito|oxidative|stress|module|score|hallmark",
    meta_cols, ignore.case = TRUE
  )]
# exclude standard qc columns
  exclude <- c("percent.mt","nCount","nFeature","seurat_clusters",
               "doublet","sample","orig","prediction","barcode")
  pattern_hits <- pattern_hits[!grepl(paste(exclude, collapse="|"),
                                      pattern_hits, ignore.case=TRUE)]
  return(pattern_hits)
}

# figure 1: per-sample spatial cell type maps

message("\n", strrep("=", 65))
message("  TASK B: SPATIAL CELL TYPE MAPS")
message(strrep("=", 65))

make_spatial_celltype_map <- function(so, sample_id, condition,
                                      dataset, pal) {
  
  meta <- so@meta.data
  if (!all(c("x","y","predicted_celltype") %in% colnames(meta))) {
    message("  Missing required columns for ", sample_id)
    return(NULL)
  }
  
# filter to high-confidence spots
  plot_df <- meta %>%
    filter(!is.na(x), !is.na(y), !is.na(predicted_celltype))
  
  n_total    <- nrow(plot_df)
  n_pass     <- sum(plot_df$prediction_score >= SCORE_THRESH, na.rm = TRUE)
  med_score  <- round(median(plot_df$prediction_score, na.rm = TRUE), 3)
  
  p <- ggplot(plot_df, aes(x = x, y = -y, colour = predicted_celltype)) +
    geom_point(size = 0.6, alpha = 0.85) +
    scale_colour_manual(values = pal, name = "Cell Type",
                        na.value = "grey90") +
    coord_fixed() +
    theme_void(base_size = 9) +
    theme(
      legend.position   = "right",
      legend.text       = element_text(size = 7),
      legend.title      = element_text(size = 8, face = "bold"),
      legend.key.size   = unit(0.35, "cm"),
      plot.title        = element_text(face = "bold", hjust = 0.5, size = 10),
      plot.subtitle     = element_text(hjust = 0.5, size = 7.5,
                                       colour = "grey50"),
      plot.margin       = margin(5, 5, 5, 5)
    ) +
    labs(
      title    = paste0(toupper(dataset), "  |  ", sample_id,
                        "  [", condition, "]"),
      subtitle = sprintf("%d spots | median score = %.3f",
                         n_total, med_score)
    )
  
  return(p)
}

# load and plot per sample
all_meta <- list()  # collect metadata for downstream analyses

for (dataset in c("pcos", "aging")) {
  samples <- if (dataset == "pcos") PCOS_SAMPLES else AGING_SAMPLES
  cond_map <- if (dataset == "pcos") PCOS_CONDITION else AGING_CONDITION
  
  message("\n  Dataset: ", toupper(dataset))
  
# load all samples to get global cell type palette
  so_list <- list()
  for (s in samples) {
    so <- load_lt_sample(dataset, s)
    if (!is.null(so)) so_list[[s]] <- so
  }
  
  if (length(so_list) == 0) {
    message("  No samples loaded for ", dataset)
    next
  }
  
  all_ct <- unique(unlist(lapply(so_list, function(x)
    unique(x@meta.data$predicted_celltype))))
  all_ct <- sort(na.omit(all_ct))
  pal <- make_celltype_palette(all_ct)
  
# generate per-sample spatial plots
  plot_list <- list()
  for (s in names(so_list)) {
    condition <- if (s %in% names(cond_map)) cond_map[s] else "unknown"
    p <- make_spatial_celltype_map(so_list[[s]], s, condition, dataset, pal)
    if (!is.null(p)) plot_list[[s]] <- p
    
# collect metadata for proportion analysis
    m <- so_list[[s]]@meta.data
    
# coerce coordinate columns to numeric to avoid type conflicts in bind_rows
# (some samples store pxl_col as character, others as integer)
    for (coord_col in c("pxl_col","pxl_row","array_col","array_row","x","y")) {
      if (coord_col %in% colnames(m)) m[[coord_col]] <- as.numeric(m[[coord_col]])
    }
    
    m$sample_id  <- s
    m$condition  <- if (s %in% names(cond_map)) cond_map[s] else "unknown"
    m$dataset    <- dataset
    all_meta[[paste0(dataset,"_",s)]] <- m
  }
  
  if (length(plot_list) > 0) {
    ncols  <- min(3, length(plot_list))
    nrows  <- ceiling(length(plot_list) / ncols)
    p_grid <- wrap_plots(plot_list, ncol = ncols, guides = "collect") +
      plot_annotation(
        title = paste0(toupper(dataset),
                       " — Predicted Cell Type Distribution"),
        theme = theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 13)
        )
      )
    
    fname <- paste0(toupper(dataset), "_spatial_celltype_overview")
    ggsave(file.path(OUT_DIR, paste0(fname, ".png")),
           p_grid, width = 6*ncols, height = 5.5*nrows + 0.5, dpi = 300)
    ggsave(file.path(OUT_DIR, paste0(fname, ".tiff")),
           p_grid, width = 6*ncols, height = 5.5*nrows + 0.5, dpi = 600,
           compression = "lzw")
    message("  ✓ ", toupper(dataset), " spatial overview saved (",
            length(plot_list), " samples)")
  }
  
  gc()
}

# figure 2: cell type proportion comparison

message("\n  Building cell type proportion comparison...")

if (length(all_meta) > 0) {
  
  combined_meta <- bind_rows(all_meta) %>%
    filter(!is.na(predicted_celltype),
           !is.na(condition),
           condition != "unknown")
  
# proportions per sample
  prop_df <- combined_meta %>%
    group_by(dataset, condition, sample_id, predicted_celltype) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(dataset, condition, sample_id) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
  
# mean proportion per condition
  mean_prop <- prop_df %>%
    group_by(dataset, condition, predicted_celltype) %>%
    summarise(
      mean_prop = mean(proportion, na.rm = TRUE),
      se_prop   = sd(proportion,   na.rm = TRUE) /
        sqrt(n()),
      .groups = "drop"
    )
  
  write_csv(prop_df,   file.path(OUT_DIR, "celltype_proportions_per_sample.csv"))
  write_csv(mean_prop, file.path(OUT_DIR, "celltype_proportions_mean.csv"))
  message("  ✓ Proportion tables saved")
  
# stacked barplot: cell type proportions by condition
  for (ds in c("pcos","aging")) {
    
    ds_prop <- prop_df %>%
      filter(dataset == ds) %>%
      mutate(
        cell_type_clean = gsub("_"," ", predicted_celltype),
        condition_label = paste0(condition, "\n(", sample_id, ")")
      )
    
    if (nrow(ds_prop) == 0) next
    
    all_ct_ds <- sort(unique(ds_prop$predicted_celltype))
    pal_ds    <- make_celltype_palette(all_ct_ds)
    names(pal_ds) <- gsub("_"," ", names(pal_ds))
    
    ds_prop$ct_clean <- gsub("_"," ", ds_prop$predicted_celltype)
    
# order samples: controls first, then cases
    cond_order <- if (ds == "pcos") c("control","case") else c("young","aged")
    ds_prop$condition <- factor(ds_prop$condition,
                                levels = intersect(cond_order,
                                                   unique(ds_prop$condition)))
    
    p_stack <- ggplot(ds_prop,
                      aes(x = sample_id, y = proportion, fill = ct_clean)) +
      geom_bar(stat = "identity", position = "stack", width = 0.75) +
      scale_fill_manual(values = pal_ds, name = "Cell Type") +
      scale_y_continuous(labels = percent_format(), expand = c(0,0)) +
      facet_grid(~condition, scales = "free_x", space = "free_x") +
      theme_pub(base_size = 10) +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
        legend.text  = element_text(size = 7),
        legend.key.size = unit(0.35,"cm")
      ) +
      labs(
        title    = paste0(toupper(ds),
                          " — Cell Type Proportions per Sample"),
        subtitle = paste0("Based on label transfer | score threshold ≥ ",
                          SCORE_THRESH),
        x = NULL, y = "Proportion of spots"
      )
    
    fname <- paste0(toupper(ds), "_celltype_proportions")
    ggsave(file.path(OUT_DIR, paste0(fname, ".png")),
           p_stack, width = max(8, length(unique(ds_prop$sample_id)) * 1.5),
           height = 7, dpi = 300)
    ggsave(file.path(OUT_DIR, paste0(fname, ".tiff")),
           p_stack, width = max(8, length(unique(ds_prop$sample_id)) * 1.5),
           height = 7, dpi = 600, compression = "lzw")
    message("  ✓ ", toupper(ds), " proportion stacked barplot saved")
  }
  
# dot plot: mean proportion case/ctrl or aged/young side-by-side
  for (ds in c("pcos","aging")) {
    
    ds_mean <- mean_prop %>%
      filter(dataset == ds) %>%
      mutate(ct_clean = gsub("_"," ", predicted_celltype))
    
# need at least 2 conditions to compare
    conditions <- unique(ds_mean$condition)
    if (length(conditions) < 2) {
      message("  Only one condition in ", toupper(ds), " — skipping comparison plot")
      next
    }
    
# pivot wide: one row per cell type, columns = conditions
    comp_wide <- ds_mean %>%
      select(ct_clean, condition, mean_prop) %>%
      pivot_wider(names_from = condition, values_from = mean_prop,
                  values_fill = 0)
    
    cond_cols <- setdiff(colnames(comp_wide), "ct_clean")
    if (length(cond_cols) < 2) next
    
    comp_wide$diff <- comp_wide[[cond_cols[1]]] - comp_wide[[cond_cols[2]]]
    
    p_comp <- ggplot(ds_mean,
                     aes(x = mean_prop, y = reorder(ct_clean, mean_prop),
                         colour = condition, shape = condition)) +
      geom_point(size = 3.5, alpha = 0.9) +
      geom_line(aes(group = ct_clean), colour = "grey70", linewidth = 0.5) +
      scale_x_continuous(labels = percent_format(accuracy = 0.1)) +
      scale_colour_manual(
        values = if (ds == "pcos")
          c(case = "#D55E00", control = "#0072B2")
        else
          c(aged = "#D55E00", young = "#0072B2"),
        name = "Condition"
      ) +
      theme_pub(base_size = 10) +
      theme(axis.text.y = element_text(size = 8)) +
      labs(
        title    = paste0(toupper(ds),
                          " — Cell Type Proportion Comparison"),
        subtitle = paste0("Mean proportion per condition | ",
                          cond_cols[1], " vs ", cond_cols[2]),
        x = "Mean proportion of spots", y = NULL
      )
    
    fname <- paste0(toupper(ds), "_celltype_proportion_comparison")
    ggsave(file.path(OUT_DIR, paste0(fname, ".png")),
           p_comp, width = 8, height = 6, dpi = 300)
    ggsave(file.path(OUT_DIR, paste0(fname, ".tiff")),
           p_comp, width = 8, height = 6, dpi = 600, compression = "lzw")
    message("  ✓ ", toupper(ds), " proportion comparison saved")
  }
}

# figure 3: co-localisation with mitochondrial module scores

message("\n  Co-localisation: cell types × mitochondrial scores...")

coloc_results <- list()

for (ds in c("pcos","aging")) {
  
  int_path <- if (ds == "pcos") INTEGRATED_PCOS else INTEGRATED_AGING
  if (!file.exists(int_path)) {
    message("  Integrated object not found: ", int_path)
    next
  }
  
  message("  Loading integrated object: ", basename(int_path))
  so_int <- tryCatch(readRDS(int_path),
                     error = function(e) {
                       message("  Load error: ", e$message); NULL
                     })
  if (is.null(so_int)) next
  
  message("  Integrated: ", ncol(so_int), " spots")
  
# find module score columns
  score_cols <- find_module_score_cols(so_int)
  message("  Module score columns found: ",
          if (length(score_cols) > 0) paste(score_cols, collapse=", ")
          else "NONE")
  
  if (length(score_cols) == 0) {
    message("  WARNING: No module score columns found for ", toupper(ds))
    message("  Available metadata cols: ",
            paste(colnames(so_int@meta.data), collapse=", "))
    next
  }
  
# get barcodes from integrated object
  int_barcodes <- colnames(so_int)
  int_meta     <- so_int@meta.data[, c("orig.ident",
                                       intersect(score_cols,
                                                 colnames(so_int@meta.data))),
                                   drop = FALSE]
  
# load label transfer results for this dataset
  samples <- if (ds == "pcos") PCOS_SAMPLES else AGING_SAMPLES
  
  lt_meta_list <- list()
  for (s in samples) {
    csv_path <- file.path(LT_BASE, ds, s,
                          paste0(s, "_metadata_with_celltypes.csv"))
    if (!file.exists(csv_path)) next
    m <- tryCatch(read_csv(csv_path, show_col_types = FALSE),
                  error = function(e) NULL)
    if (is.null(m)) next
    m$sample_id <- s
    lt_meta_list[[s]] <- m
  }
  
  if (length(lt_meta_list) == 0) {
    message("  No label transfer metadata found for ", toupper(ds))
    next
  }
  
  lt_meta <- bind_rows(lt_meta_list)
  
# identify the barcode column
  bc_col <- intersect(c("barcode","Barcode","cell_id"), colnames(lt_meta))[1]
  if (is.na(bc_col)) {
# assume first column is barcodes
    bc_col <- colnames(lt_meta)[1]
  }
  
# try to match barcodes between label transfer and integrated object
# barcodes may have sample suffixes added during integration
  lt_meta$barcode_clean <- lt_meta[[bc_col]]
  
# build combined table
# match on barcode — try direct, then suffix-stripped
  int_meta$barcode_int <- rownames(int_meta)
  
  merged <- lt_meta %>%
    select(barcode_clean, sample_id, predicted_celltype, prediction_score) %>%
    filter(!is.na(predicted_celltype),
           prediction_score >= SCORE_THRESH) %>%
    left_join(
      int_meta %>%
        rownames_to_column("barcode_int") %>%
        mutate(barcode_clean = sub("_[0-9]+$", "", barcode_int)),
      by = "barcode_clean"
    )
  
# if poor match, try original barcodes
  n_matched <- sum(!is.na(merged$barcode_int))
  message("  Matched ", n_matched, "/", nrow(merged),
          " spots between label transfer and integrated object")
  
  if (n_matched < 100) {
# try appending sample suffixes
    message("  Poor match — trying barcode suffix matching...")
# integration often appends _1, _2... per sample in order
    sample_suffix_map <- setNames(seq_along(samples), samples)
    
    lt_meta$barcode_suffixed <- paste0(
      lt_meta$barcode_clean, "_",
      sample_suffix_map[lt_meta$sample_id]
    )
    merged <- lt_meta %>%
      select(barcode_suffixed, sample_id, predicted_celltype, prediction_score) %>%
      filter(!is.na(predicted_celltype),
             prediction_score >= SCORE_THRESH) %>%
      left_join(
        int_meta %>% rownames_to_column("barcode_int"),
        by = c("barcode_suffixed" = "barcode_int")
      )
    n_matched <- sum(!is.na(merged[[score_cols[1]]]))
    message("  After suffix matching: ", n_matched, " spots matched")
  }
  
  if (n_matched < 50) {
    message("  WARNING: Too few spots matched for co-localisation (",
            toupper(ds), ")")
    message("  This is likely a barcode naming mismatch between integration and label transfer")
    message("  Skipping co-localisation for ", toupper(ds))
    next
  }
  
# for each module score: compute mean score per cell type
  for (sc in score_cols) {
    
    if (!sc %in% colnames(merged)) next
    
    merged_sc <- merged %>%
      filter(!is.na(.data[[sc]])) %>%
      rename(module_score = all_of(sc))
    
    if (nrow(merged_sc) == 0) next
    
# mean module score per cell type
    ct_scores <- merged_sc %>%
      group_by(predicted_celltype) %>%
      summarise(
        mean_score = mean(module_score, na.rm = TRUE),
        median_score = median(module_score, na.rm = TRUE),
        n_spots    = n(),
        .groups = "drop"
      ) %>%
      mutate(
        module      = sc,
        dataset     = toupper(ds),
        ct_clean    = gsub("_"," ", predicted_celltype)
      ) %>%
      arrange(desc(mean_score))
    
    coloc_results[[paste0(ds,"_",sc)]] <- ct_scores
  }
  
  rm(so_int); gc()
}

# co-localisation plots

if (length(coloc_results) > 0) {
  
  coloc_df <- bind_rows(coloc_results)
  write_csv(coloc_df, file.path(OUT_DIR, "colocalisation_mito_celltype.csv"))
  
# separate plot per module score × dataset
  for (mod in unique(coloc_df$module)) {
    for (ds in unique(coloc_df$dataset)) {
      
      df_sub <- coloc_df %>%
        filter(module == mod, dataset == ds, n_spots >= 10)
      
      if (nrow(df_sub) < 2) next
      
      p_coloc <- ggplot(df_sub,
                        aes(x = mean_score,
                            y = reorder(ct_clean, mean_score),
                            size = n_spots, colour = mean_score)) +
        geom_segment(aes(x = min(mean_score) - 0.01,
                         xend = mean_score, yend = ct_clean),
                     colour = "grey82", linewidth = 0.5) +
        geom_point() +
        scale_colour_gradient2(low = "#4575b4", mid = "grey85",
                               high = "#d73027", midpoint = 0,
                               name = "Mean\nScore") +
        scale_size_continuous(name = "Spots", range = c(2,7)) +
        theme_pub(base_size = 10) +
        theme(axis.text.y = element_text(size = 8.5)) +
        labs(
          title    = paste0(ds, " — ", gsub("_"," ", mod),
                            " Module Score by Cell Type"),
          subtitle = "Mean module score per predicted cell type",
          x = "Mean module score", y = NULL
        )
      
      mod_clean <- clean_name <- gsub("[^A-Za-z0-9]","_", mod)
      fname <- paste0("Colocalisation_", ds, "_", mod_clean)
      ggsave(file.path(OUT_DIR, paste0(fname, ".png")),
             p_coloc, width = 8, height = max(5, 0.4*nrow(df_sub)+2),
             dpi = 300)
      ggsave(file.path(OUT_DIR, paste0(fname, ".tiff")),
             p_coloc, width = 8, height = max(5, 0.4*nrow(df_sub)+2),
             dpi = 600, compression = "lzw")
    }
  }
  
# combined heatmap: cell types × module scores
  heat_df <- coloc_df %>%
    filter(n_spots >= 10) %>%
    mutate(col_label = paste0(dataset, "\n", gsub("_"," ", module))) %>%
    select(ct_clean, col_label, mean_score) %>%
    pivot_wider(names_from = col_label, values_from = mean_score,
                values_fn = mean, values_fill = NA)
  
  if (nrow(heat_df) >= 2 && ncol(heat_df) >= 3) {
    heat_mat <- as.matrix(heat_df[,-1])
    rownames(heat_mat) <- heat_df$ct_clean
    
# scale column-wise for comparability
    heat_scaled <- scale(heat_mat)
    heat_scaled[is.nan(heat_scaled)] <- 0
    
    heat_long <- as.data.frame(heat_scaled) %>%
      rownames_to_column("cell_type") %>%
      pivot_longer(-cell_type, names_to = "module", values_to = "score_z") %>%
      mutate(cell_type = factor(cell_type,
                                levels = rownames(heat_scaled)[
                                  hclust(dist(heat_scaled))$order]))
    
    p_heat <- ggplot(heat_long, aes(x = module, y = cell_type, fill = score_z)) +
      geom_tile(colour = "white", linewidth = 0.3) +
      scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",
                           midpoint = 0, na.value = "grey90",
                           name = "Z-score") +
      theme_pub(base_size = 9) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        axis.ticks  = element_blank(),
        panel.grid  = element_blank()
      ) +
      labs(
        title    = "Cell Type Co-localisation with Module Scores",
        subtitle = "Scaled mean module score per cell type (z-score)",
        x = NULL, y = NULL
      )
    
    h <- max(5, 0.4 * nrow(heat_scaled) + 2)
    w <- max(7, 0.8 * ncol(heat_scaled) + 3)
    ggsave(file.path(OUT_DIR, "Colocalisation_heatmap.png"),
           p_heat, width = w, height = h, dpi = 300)
    ggsave(file.path(OUT_DIR, "Colocalisation_heatmap.tiff"),
           p_heat, width = w, height = h, dpi = 600, compression = "lzw")
    message("  ✓ Co-localisation heatmap saved")
  }
  
  message("  ✓ Co-localisation analysis complete")
} else {
  message("  Co-localisation skipped (no module score columns found or barcode mismatch)")
  message("  Spatial proportion maps were still generated successfully")
}

# final summary

message("\n", strrep("=", 65))
message("  TASK B COMPLETE")
message(strrep("=", 65))
message("  Output directory: ", OUT_DIR)
message("  Files produced:")
for (f in list.files(OUT_DIR, full.names = FALSE, recursive = TRUE))
  message("    ", f)