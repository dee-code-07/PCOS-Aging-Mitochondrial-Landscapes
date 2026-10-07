#!/usr/bin/env Rscript
# =============================================================================
# Script: 28_slingshot_multilineage.R
# Purpose: Document ALL Slingshot lineages before selection (Minor Concern 5)
#          - Diagnostic: verify lineage curve endpoints match stated paths
#          - Checkpoint: save gran/sce/sds for cheap future reruns
#          - Output: publication-quality PDF (rasterized UMAP, vector curves)
# Output: /home/deekshah/mini_project/analysis/28_slingshot_multilineage/
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(slingshot)
  library(SingleCellExperiment)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
  library(ggrastr)
})

set.seed(42)
options(future.globals.maxSize = Inf)

# --- Configuration ---
base_dir <- "/home/deekshah/mini_project"
out_dir  <- file.path(base_dir, "analysis/28_slingshot_multilineage")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pcos_path  <- file.path(base_dir, "scrna/output/07_annotation/pcos/pcos_annotated.rds")
aging_path <- file.path(base_dir, "scrna/output/07_annotation/aging/aging_annotated.rds")

gran_types <- c("granulosa", "granulosa_antral", "granulosa_luteal", "granulosa_preantral")

# Open a diagnostic log file (append mode — one file for both datasets)
diag_file <- file.path(out_dir, "diagnostic_output.txt")
diag_con  <- file(diag_file, open = "wt")

diag_log <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  writeLines(msg, diag_con)
}

# --- Helper: process one dataset ---
process_dataset <- function(obj_path, dataset_name, out_dir) {

  diag_log("")
  diag_log(sprintf("========== Processing %s ==========", dataset_name))
  obj <- readRDS(obj_path)

  # Add condition
  if (dataset_name == "PCOS") {
    obj$condition <- ifelse(grepl("case", obj$sample_id, ignore.case = TRUE), "PCOS",
                     ifelse(grepl("ctrl|control", obj$sample_id, ignore.case = TRUE), "Control", NA))
  } else {
    obj$condition <- ifelse(grepl("_a_|aged", obj$sample_id, ignore.case = TRUE), "Aged",
                     ifelse(grepl("_y_|young", obj$sample_id, ignore.case = TRUE), "Young", NA))
  }

  # Subset to granulosa
  avail_types <- intersect(gran_types, unique(obj$final_celltype))
  diag_log("Available granulosa types: ", paste(avail_types, collapse = ", "))

  if (length(avail_types) == 0) {
    diag_log("No granulosa cells found, skipping.")
    return(NULL)
  }

  gran <- subset(obj, final_celltype %in% avail_types)
  total_gran <- ncol(gran)
  rm(obj); gc()

  diag_log("Total granulosa cells: ", total_gran)

  # Re-process
  DefaultAssay(gran) <- "RNA"
  tryCatch(gran <- JoinLayers(gran), error = function(e) NULL)
  gran <- NormalizeData(gran, verbose = FALSE)
  gran <- FindVariableFeatures(gran, nfeatures = 3000, verbose = FALSE)
  gran <- ScaleData(gran, verbose = FALSE)
  gran <- RunPCA(gran, npcs = 30, verbose = FALSE)
  gran <- RunUMAP(gran, dims = 1:20, verbose = FALSE)

  # Convert to SCE
  sce <- as.SingleCellExperiment(gran)
  reducedDim(sce, "UMAP") <- Embeddings(gran, "umap")

  # Determine start cluster
  start_clus <- if ("granulosa_preantral" %in% avail_types) "granulosa_preantral" else avail_types[1]
  diag_log("Start cluster: ", start_clus)

  # Run Slingshot
  sce <- slingshot(sce, clusterLabels = "final_celltype", reducedDim = "UMAP",
                   start.clus = start_clus)

  # Extract lineage info
  sds <- SlingshotDataSet(sce)
  lineages <- slingLineages(sds)
  n_lineages <- length(lineages)
  diag_log("Number of lineages: ", n_lineages)

  # =========================================================================
  # CHECKPOINT: save gran, sce, sds for cheap future reruns (Step 3)
  # =========================================================================
  checkpoint_path <- file.path(out_dir, paste0(dataset_name, "_checkpoint.rds"))
  diag_log("Saving checkpoint to: ", checkpoint_path)
  saveRDS(list(gran = gran, sce = sce, sds = sds),
          file = checkpoint_path)
  diag_log("Checkpoint saved.")

  # =========================================================================
  # DIAGNOSTICS (Step 1): verify lineage curves match stated paths
  # =========================================================================
  diag_log("")
  diag_log("--- DIAGNOSTICS ---")

  # 1. Lineage cluster order
  diag_log("1. Lineage cluster order (slingLineages):")
  for (i in seq_along(lineages)) {
    diag_log("   Lineage ", i, ": ", paste(lineages[[i]], collapse = " -> "))
  }

  # 2. Curve count and names
  curves <- slingCurves(sds)
  diag_log("2. Number of curves: ", length(curves))
  diag_log("   Curve names: ", paste(names(curves), collapse = ", "))

  # 3. For each lineage curve, find which cell type cluster its ENDPOINT lands in
  diag_log("3. Curve endpoint nearest cell types (top 20 nearest cells):")
  umap_coords <- Embeddings(gran, "umap")

  for (i in seq_along(curves)) {
    curve_data <- curves[[i]]
    end_point <- curve_data$s[tail(curve_data$ord, 1), ]
    dists <- sqrt((umap_coords[, 1] - end_point[1])^2 +
                  (umap_coords[, 2] - end_point[2])^2)
    nearest_cells <- order(dists)[1:20]
    ct_table <- table(gran$final_celltype[nearest_cells])

    diag_log("   Lineage ", i, " curve endpoint (UMAP: ",
             round(end_point[1], 2), ", ", round(end_point[2], 2), "):")
    for (ct in names(ct_table)) {
      diag_log("     ", ct, ": ", ct_table[ct])
    }

    # Also report the START point nearest cell types
    start_point <- curve_data$s[curve_data$ord[1], ]
    dists_start <- sqrt((umap_coords[, 1] - start_point[1])^2 +
                        (umap_coords[, 2] - start_point[2])^2)
    nearest_start <- order(dists_start)[1:20]
    ct_start <- table(gran$final_celltype[nearest_start])
    diag_log("   Lineage ", i, " curve start (UMAP: ",
             round(start_point[1], 2), ", ", round(start_point[2], 2), "):")
    for (ct in names(ct_start)) {
      diag_log("     ", ct, ": ", ct_start[ct])
    }
  }

  # 4. Cross-check pseudotime-based cell type composition for each lineage
  diag_log("4. Cell type composition per pseudotime assignment:")
  for (i in seq_len(n_lineages)) {
    pt_col <- paste0("slingPseudotime_", i)
    if (pt_col %in% colnames(colData(sce))) {
      pt_cells <- which(!is.na(colData(sce)[[pt_col]]))
      ct_table <- table(gran$final_celltype[pt_cells])
      diag_log("   ", pt_col, " (n = ", length(pt_cells), " cells):")
      for (ct in names(ct_table)) {
        diag_log("     ", ct, ": ", ct_table[ct],
                 " (", round(100 * ct_table[ct] / length(pt_cells), 1), "%)")
      }

      # Condition composition for this lineage
      cond_table <- table(gran$condition[pt_cells])
      diag_log("   Condition breakdown:")
      for (cond in names(cond_table)) {
        diag_log("     ", cond, ": ", cond_table[cond],
                 " (", round(100 * cond_table[cond] / length(pt_cells), 1), "%)")
      }
    }
  }

  # 5. Exclusion statistics: cells NOT in the primary lineage
  diag_log("5. Exclusion analysis (cells not in Lineage 1):")
  pt1_cells <- which(!is.na(colData(sce)$slingPseudotime_1))
  n_primary <- length(pt1_cells)
  n_excluded <- total_gran - n_primary
  diag_log("   Total granulosa cells: ", total_gran)
  diag_log("   Cells in primary lineage (Lineage 1): ", n_primary,
           " (", round(100 * n_primary / total_gran, 1), "%)")
  diag_log("   Cells excluded by retaining only primary: ", n_excluded,
           " (", round(100 * n_excluded / total_gran, 1), "%)")

  # Condition breakdown of excluded cells
  if (n_excluded > 0) {
    # Cells in any lineage vs not in lineage 1
    all_pt_cols <- paste0("slingPseudotime_", seq_len(n_lineages))
    in_any <- rep(FALSE, total_gran)
    for (col in all_pt_cols) {
      if (col %in% colnames(colData(sce))) {
        in_any <- in_any | !is.na(colData(sce)[[col]])
      }
    }
    in_primary <- !is.na(colData(sce)$slingPseudotime_1)

    # Cells only in non-primary lineages (in a secondary lineage but NOT in lineage 1)
    only_secondary <- in_any & !in_primary
    n_only_secondary <- sum(only_secondary)
    diag_log("   Cells exclusively in non-primary lineage(s): ", n_only_secondary)
    if (n_only_secondary > 0) {
      cond_excl <- table(gran$condition[only_secondary])
      diag_log("   Condition breakdown of exclusively non-primary cells:")
      for (cond in names(cond_excl)) {
        diag_log("     ", cond, ": ", cond_excl[cond],
                 " (", round(100 * cond_excl[cond] / n_only_secondary, 1), "%)")
      }
    }

    # Cells in NO lineage at all
    in_none <- !in_any
    n_none <- sum(in_none)
    diag_log("   Cells not assigned to ANY lineage: ", n_none,
             " (", round(100 * n_none / total_gran, 1), "%)")
    if (n_none > 0) {
      cond_none <- table(gran$condition[in_none])
      ct_none <- table(gran$final_celltype[in_none])
      diag_log("   Cell types of unassigned cells:")
      for (ct in names(ct_none)) {
        diag_log("     ", ct, ": ", ct_none[ct])
      }
      diag_log("   Condition of unassigned cells:")
      for (cond in names(cond_none)) {
        diag_log("     ", cond, ": ", cond_none[cond])
      }
    }
  }

  # =========================================================================
  # LINEAGE STATISTICS TABLE (for CSV output)
  # =========================================================================
  lineage_stats <- data.frame()

  for (i in seq_len(n_lineages)) {
    lin_name <- paste0("Lineage", i)
    pt_col <- paste0("slingPseudotime_", i)

    if (pt_col %in% colnames(colData(sce))) {
      pt_vals <- colData(sce)[[pt_col]]
      has_pt <- !is.na(pt_vals)
      n_cells <- sum(has_pt)

      # Condition composition
      cond_table <- table(gran$condition[has_pt])
      cond_props <- prop.table(cond_table)

      # Cell type composition
      ct_table <- table(gran$final_celltype[has_pt])

      for (cond in names(cond_table)) {
        lineage_stats <- rbind(lineage_stats, data.frame(
          dataset = dataset_name,
          lineage = lin_name,
          path = paste(lineages[[i]], collapse = " -> "),
          total_cells = n_cells,
          total_granulosa = total_gran,
          pct_of_granulosa = round(100 * n_cells / total_gran, 1),
          condition = cond,
          n_condition = as.integer(cond_table[cond]),
          prop_condition = round(as.numeric(cond_props[cond]), 4),
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Save lineage stats
  write_csv(lineage_stats,
            file.path(out_dir, paste0(dataset_name, "_lineage_statistics.csv")))

  # =========================================================================
  # FIGURES: Publication-quality PDF (rasterized UMAP, vector curves/text)
  # =========================================================================
  umap_df <- data.frame(
    UMAP_1 = Embeddings(gran, "umap")[, 1],
    UMAP_2 = Embeddings(gran, "umap")[, 2],
    celltype = gran$final_celltype,
    condition = gran$condition
  )

  # --- Cell-type color palette (consistent, colorblind-friendly) ---
  celltype_colors <- c(
    "granulosa"           = "#F28E2B",
    "granulosa_antral"    = "#4E79A7",
    "granulosa_luteal"    = "#59A14F",
    "granulosa_preantral" = "#B07AA1"
  )
  # Keep only colors for types actually present
  celltype_colors <- celltype_colors[names(celltype_colors) %in% unique(umap_df$celltype)]

  # --- Panel A: UMAP colored by cell type with all lineage curves ---
  p_celltype <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
    rasterise(geom_point(size = 0.3, alpha = 0.4), dpi = 300) +
    scale_color_manual(values = celltype_colors,
                       labels = gsub("_", " ", names(celltype_colors))) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey40"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    ) +
    labs(title = paste0(dataset_name, " — All Slingshot lineages"),
         subtitle = paste0(n_lineages, " lineages inferred from ",
                           format(total_gran, big.mark = ","), " granulosa cells"),
         color = "Cell type") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1),
                                nrow = 1))

  # Add curves for each lineage
  curve_colors <- c("#E63946", "#457B9D", "#2A9D8F", "#E9C46A", "#F4A261")

  for (i in seq_len(n_lineages)) {
    curve_data <- slingCurves(sds)[[i]]
    curve_df <- as.data.frame(curve_data$s[curve_data$ord, ])
    colnames(curve_df) <- c("UMAP_1", "UMAP_2")

    # Endpoint for label placement (with slight offset to avoid overlap)
    end_x <- tail(curve_df$UMAP_1, 1)
    end_y <- tail(curve_df$UMAP_2, 1)

    p_celltype <- p_celltype +
      geom_path(data = curve_df, aes(x = UMAP_1, y = UMAP_2),
                color = curve_colors[min(i, length(curve_colors))],
                linewidth = 1.2, inherit.aes = FALSE) +
      annotate("label", x = end_x, y = end_y,
               label = paste0("L", i), fontface = "bold",
               color = curve_colors[min(i, length(curve_colors))],
               fill = "white", size = 4, label.size = 0.3,
               label.padding = unit(0.2, "lines"))
  }

  # --- Panel B: UMAP colored by condition ---
  if (dataset_name == "PCOS") {
    cond_colors <- c("Control" = "#457B9D", "PCOS" = "#E63946")
  } else {
    cond_colors <- c("Young" = "#457B9D", "Aged" = "#E63946")
  }

  p_condition <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = condition)) +
    rasterise(geom_point(size = 0.3, alpha = 0.4), dpi = 300) +
    scale_color_manual(values = cond_colors) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    ) +
    labs(title = paste0(dataset_name, " — Condition overlay"),
         color = "Condition") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

  # --- Panel C: Pseudotime density per lineage split by condition ---
  pt_long <- data.frame()
  for (i in seq_len(n_lineages)) {
    pt_col <- paste0("slingPseudotime_", i)
    if (pt_col %in% colnames(colData(sce))) {
      pt_vals <- colData(sce)[[pt_col]]
      valid <- !is.na(pt_vals)
      pt_long <- rbind(pt_long, data.frame(
        pseudotime = pt_vals[valid],
        lineage = paste0("Lineage ", i, "\n(",
                         paste(lineages[[i]], collapse = " → "), ")"),
        condition = gran$condition[valid]
      ))
    }
  }

  if (dataset_name == "PCOS") {
    fill_colors <- c("Control" = "#457B9D", "PCOS" = "#E63946")
  } else {
    fill_colors <- c("Young" = "#457B9D", "Aged" = "#E63946")
  }

  p_density <- ggplot(pt_long, aes(x = pseudotime, fill = condition)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ lineage, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = fill_colors) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank()
    ) +
    labs(title = "Pseudotime density by lineage and condition",
         x = "Pseudotime", y = "Density", fill = "Condition")

  # --- Combine panels ---
  combined <- (p_celltype | p_condition) / p_density +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(
      title = paste0("Supplementary Figure: ", dataset_name,
                     " Slingshot lineage analysis"),
      theme = theme(plot.title = element_text(face = "bold", size = 15))
    )

  # --- Save as PDF (cairo_pdf for proper alpha rendering) ---
  ggsave(file.path(out_dir, paste0(dataset_name, "_all_lineages.pdf")),
         combined, width = 16, height = 14, device = cairo_pdf)

  diag_log("")
  diag_log(sprintf("Done with %s. Outputs saved to %s", dataset_name, out_dir))

  rm(gran, sce, sds); gc()

  return(lineage_stats)
}

# --- Main ---
diag_log("Starting Slingshot multilineage analysis")
diag_log(paste0("Timestamp: ", Sys.time()))
diag_log("=========================================")

pcos_stats  <- process_dataset(pcos_path, "PCOS", out_dir)
aging_stats <- process_dataset(aging_path, "Aging", out_dir)

# Combine and save
all_stats <- bind_rows(pcos_stats, aging_stats)
write_csv(all_stats, file.path(out_dir, "all_lineage_statistics_combined.csv"))

diag_log("")
diag_log("========== COMPLETE ==========")
diag_log(paste0("Timestamp: ", Sys.time()))
diag_log("All outputs saved to: ", out_dir)
diag_log("Files generated:")
diag_log("  - PCOS_checkpoint.rds")
diag_log("  - Aging_checkpoint.rds")
diag_log("  - PCOS_lineage_statistics.csv")
diag_log("  - Aging_lineage_statistics.csv")
diag_log("  - all_lineage_statistics_combined.csv")
diag_log("  - PCOS_all_lineages.pdf")
diag_log("  - Aging_all_lineages.pdf")
diag_log("  - diagnostic_output.txt")

close(diag_con)
