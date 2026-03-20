# script 02: hdwgcna core tutorial plots
# project: single-cell and spatial landscapes of mitochondrial dysfunction
# in pcos and ovarian aging
# purpose:
# generates all plots specified in the official hdwgcna tutorial:
# https://smorabit.github.io/hdwgcna/articles/basic_tutorial.html
# run after the main hdwgcna script (01_hdwgcna_pcos_aging.r) completes.
# plots generated per cell type:
# p1  - plotdendrogram           (gene dendrogram + module color bar)
# p2  - plotkmes                 (genes ranked by kme per module)
# p3  - modulefeatureplot hmes   (harmonized me on full umap)
# p4  - modulefeatureplot scores (ucell hub gene scores on umap)
# p5  - moduleradarplot          (module activity across conditions)
# p6  - modulecorrelogram        (module-module correlation matrix)
# p7  - dotplot hmes             (module eigengenes × condition)
# p8  - gethubgenes table        (saved as csv)
# p9  - moduleexprscore umap     (hub gene scores ucell)
# p10 - getmodules table         (full module assignment saved)
# references:
# morabito s et al. cell rep methods. 2023.
# https://doi.org/10.1016/j.crmeth.2023.100498
# langfelder p & horvath s. bmc bioinformatics. 2008.
# https://doi.org/10.1186/1471-2105-9-559

# libraries & paths
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(ggplot2)
  library(viridis)
  library(RColorBrewer)
  library(scales)
  library(WGCNA)
  library(hdWGCNA)
  library(UCell)
  library(corrplot)   # for ModuleCorrelogram
})

# install corrplot if missing
if (!requireNamespace("corrplot", quietly = TRUE)) {
  install.packages("corrplot")
  library(corrplot)
}
if (!requireNamespace("ggradar", quietly = TRUE)) {
  devtools::install_github("ricardo-bion/ggradar")
}

theme_set(theme_cowplot(font_size = 12))
set.seed(42)

# paths
base_dir  <- "E:/Documents/mini_project"
out_dir   <- file.path(base_dir, "hdwgcna")
plot_dir  <- file.path(out_dir, "plots")
rds_dir   <- file.path(out_dir, "rds")
tut_dir   <- file.path(plot_dir, "tutorial_plots")  # dedicated folder

# create tutorial plots subfolder per cell type
for (ct in c("granulosa", "stromal", "immune")) {
  d <- file.path(tut_dir, ct)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# colors (consistent with script 01)
condition_colors <- c(
  "PCOS_Case"    = "#E63946",
  "PCOS_Control" = "#F4A261",
  "Aged"         = "#457B9D",
  "YoungControl" = "#A8DADC"
)

broad_colors <- c(
  "Granulosa" = "#E07B54",
  "Stromal"   = "#5B8DB8",
  "Immune"    = "#6AB187",
  "Oocyte"    = "#C49BBB",
  "Other"     = "#AAAAAA"
)

# load completed hdwgcna object
cat("Loading completed hdWGCNA Seurat object...\n")

# load the most complete checkpoint available
rds_candidates <- c(
  file.path(rds_dir, "seurat_hdwgcna_all_complete.rds"),
  file.path(rds_dir, "seurat_hdwgcna_stromal.rds"),
  file.path(rds_dir, "seurat_hdwgcna_granulosa.rds")
)

rds_file <- rds_candidates[file.exists(rds_candidates)][1]

if (is.na(rds_file)) {
  stop("No hdWGCNA RDS checkpoint found in: ", rds_dir,
       "\nPlease run 01_hdWGCNA_PCOS_Aging.R first.")
}

cat("  Loading:", rds_file, "\n")
seurat_obj <- readRDS(rds_file)
cat("✓ Object loaded.\n")
cat("  Cells:", ncol(seurat_obj), "\n")
cat("  Metadata columns:", paste(colnames(seurat_obj@meta.data), collapse=", "), "\n\n")

# detect which wgcna experiments are present
available_wgcna <- names(seurat_obj@misc)[
  sapply(names(seurat_obj@misc), function(n) {
    tryCatch({
      !is.null(seurat_obj@misc[[n]]$wgcna_params)
    }, error = function(e) FALSE)
  })
]

# fallback: check common names
expected_wgcna <- c("Granulosa_WGCNA", "Stromal_WGCNA", "Immune_WGCNA")
present_wgcna  <- expected_wgcna[expected_wgcna %in% names(seurat_obj@misc)]

if (length(present_wgcna) == 0) {
# try detecting any misc slot with modules
  present_wgcna <- names(seurat_obj@misc)[
    sapply(names(seurat_obj@misc), function(n) {
      tryCatch(!is.null(seurat_obj@misc[[n]]$wgcna_modules), error=function(e) FALSE)
    })
  ]
}

cat("✓ hdWGCNA experiments found:", paste(present_wgcna, collapse=", "), "\n\n")

# master tutorial plot function
# generates every plot from the official tutorial for one cell type.

generate_tutorial_plots <- function(seurat_obj,
                                    wgcna_name,      # e.g. "Granulosa_WGCNA"
                                    cell_type_name,  # e.g. "Granulosa"
                                    ct_short,        # e.g. "granulosa"
                                    condition_col = "condition",
                                    broad_col     = "broad_celltype") {
  
  pdir <- file.path(tut_dir, ct_short)
  cat("\n", strrep("═", 65), "\n")
  cat("  Tutorial Plots:", cell_type_name, "| Experiment:", wgcna_name, "\n")
  cat(strrep("═", 65), "\n\n")
  
# verify experiment exists
  if (!wgcna_name %in% names(seurat_obj@misc)) {
    cat("  ✗ Experiment '", wgcna_name, "' not found. Skipping.\n")
    return(invisible(seurat_obj))
  }
  
# get module table
  modules <- tryCatch(
    GetModules(seurat_obj, wgcna_name = wgcna_name),
    error = function(e) {
      cat("  ✗ GetModules failed:", conditionMessage(e), "\n"); return(NULL)
    }
  )
  if (is.null(modules)) return(invisible(seurat_obj))
  
  mods     <- levels(modules$module)
  mods     <- mods[mods != "grey"]
  n_mods   <- length(mods)
  cat("  Modules found:", n_mods, "—", paste(mods, collapse=", "), "\n\n")
  
# plot p1: dendrogram
# the canonical hdwgcna dendrogram: each leaf = one gene,
# colored band at bottom = module assignment.
# grey = unassigned genes (to be ignored in interpretation).
  cat("[P1] PlotDendrogram...\n")
  tryCatch({
    pdf(file.path(pdir, paste0("P1_dendrogram_", ct_short, ".pdf")),
        width = 14, height = 6)
    PlotDendrogram(
      seurat_obj,
      main       = paste0(cell_type_name,
                          " – Gene Co-expression Dendrogram\n",
                          "(Each leaf = gene | Color band = module assignment | Grey = unassigned)"),
      wgcna_name = wgcna_name
    )
    dev.off()
    
    png(file.path(pdir, paste0("P1_dendrogram_", ct_short, ".png")),
        width = 14, height = 6, units = "in", res = 300)
    PlotDendrogram(
      seurat_obj,
      main       = paste0(cell_type_name,
                          " – Gene Co-expression Dendrogram\n",
                          "(Each leaf = gene | Color band = module assignment | Grey = unassigned)"),
      wgcna_name = wgcna_name
    )
    dev.off()
    cat("  ✓ Dendrogram saved.\n")
  }, error = function(e) {
    cat("  ✗ Dendrogram error:", conditionMessage(e), "\n")
    try(dev.off(), silent = TRUE)
  })
  
# plot p2: kme distributions per module
# plotkmes shows the ranked kme scores for each module.
# hub genes (high kme) are biologically the most representative.
# this verifies that each module has a clear "hub" structure.
  cat("[P2] PlotKMEs (kME distributions)...\n")
  tryCatch({
    ncols_kme <- min(4, n_mods)
    p_kme <- PlotKMEs(
      seurat_obj,
      ncol       = ncols_kme,
      wgcna_name = wgcna_name
    )
# add annotation
    p_kme_final <- p_kme +
      plot_annotation(
        title    = paste0(cell_type_name, " – Gene kME Distributions Per Module"),
        subtitle = paste0("kME = eigengene-based connectivity (Pearson correlation with module eigengene)\n",
                          "Hub genes have highest kME values | n = ", n_mods, " modules"),
        theme    = theme(
          plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5)
        )
      )
    ggsave(file.path(pdir, paste0("P2_kME_distributions_", ct_short, ".pdf")),
           p_kme_final,
           width  = min(20, ncols_kme * 4 + 2),
           height = ceiling(n_mods / ncols_kme) * 3 + 2,
           dpi    = 300)
    ggsave(file.path(pdir, paste0("P2_kME_distributions_", ct_short, ".png")),
           p_kme_final,
           width  = min(20, ncols_kme * 4 + 2),
           height = ceiling(n_mods / ncols_kme) * 3 + 2,
           dpi    = 300)
    cat("  ✓ kME distributions saved.\n")
  }, error = function(e) {
    cat("  ✗ PlotKMEs error:", conditionMessage(e), "\n")
  })
  
# plot p3: modulefeatureplot — hmes on umap
# hmes (harmonized module eigengenes) projected onto the full umap.
# shows which spatial regions of the umap each module dominates.
# batch-corrected, so reflects true biology not sample artefacts.
  cat("[P3] ModuleFeaturePlot — harmonized MEs on UMAP...\n")
  tryCatch({
    plot_list_hme <- ModuleFeaturePlot(
      seurat_obj,
      features   = "hMEs",
      order      = TRUE,
      wgcna_name = wgcna_name
    )
    n_per_row <- min(4, n_mods)
    p_hme <- wrap_plots(plot_list_hme, ncol = n_per_row) +
      plot_annotation(
        title    = paste0(cell_type_name,
                          " – Harmonized Module Eigengenes on UMAP"),
        subtitle = paste0("Each panel = one co-expression module | ",
                          "Color intensity = harmonized ME score per cell\n",
                          "Harmony batch-corrected across sample_id"),
        theme    = theme(
          plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5)
        )
      )
    ggsave(file.path(pdir, paste0("P3_ModuleFeaturePlot_hMEs_", ct_short, ".pdf")),
           p_hme,
           width  = n_per_row * 4,
           height = ceiling(n_mods / n_per_row) * 4,
           dpi    = 300,
           limitsize = FALSE)
    ggsave(file.path(pdir, paste0("P3_ModuleFeaturePlot_hMEs_", ct_short, ".png")),
           p_hme,
           width  = n_per_row * 4,
           height = ceiling(n_mods / n_per_row) * 4,
           dpi    = 300,
           limitsize = FALSE)
    cat("  ✓ hME feature plots saved.\n")
  }, error = function(e) {
    cat("  ✗ ModuleFeaturePlot hMEs error:", conditionMessage(e), "\n")
  })
  
# plot p4: modulefeatureplot — hub gene ucell scores
# hub gene signature scores use the top-kme genes as a gene set.
# ucell scores are rank-based and robust to dataset size differences.
# validates that hub genes capture the same biology as the eigengene.
  cat("[P4] ModuleExprScore + ModuleFeaturePlot — UCell hub gene scores...\n")
  tryCatch({
    seurat_obj <- ModuleExprScore(
      seurat_obj,
      n_genes    = 25,
      method     = "UCell",
      wgcna_name = wgcna_name
    )
    
    plot_list_scores <- ModuleFeaturePlot(
      seurat_obj,
      features   = "scores",
      order      = "shuffle",
      ucell      = TRUE,
      wgcna_name = wgcna_name
    )
    
    n_per_row <- min(4, n_mods)
    p_scores <- wrap_plots(plot_list_scores, ncol = n_per_row) +
      plot_annotation(
        title    = paste0(cell_type_name,
                          " – Hub Gene Signature Scores on UMAP"),
        subtitle = paste0("UCell method | Top 25 hub genes (highest kME) per module\n",
                          "Rank-based scoring robust to dataset size variation"),
        theme    = theme(
          plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5)
        )
      )
    ggsave(file.path(pdir, paste0("P4_ModuleFeaturePlot_scores_", ct_short, ".pdf")),
           p_scores,
           width  = n_per_row * 4,
           height = ceiling(n_mods / n_per_row) * 4,
           dpi    = 300,
           limitsize = FALSE)
    ggsave(file.path(pdir, paste0("P4_ModuleFeaturePlot_scores_", ct_short, ".png")),
           p_scores,
           width  = n_per_row * 4,
           height = ceiling(n_mods / n_per_row) * 4,
           dpi    = 300,
           limitsize = FALSE)
    cat("  ✓ Hub gene score feature plots saved.\n")
  }, error = function(e) {
    cat("  ✗ ModuleExprScore/FeaturePlot error:", conditionMessage(e), "\n")
  })
  
# plot p5: moduleradarplot — module activity across conditions
# radar plots show the relative activity of each module per condition group.
# ideal for comparing pcos_case vs pcos_control vs aged vs youngcontrol.
# shared elevated modules across pcos and aged = convergent dysfunction.
  cat("[P5] ModuleRadarPlot — module activity by condition...\n")
  tryCatch({
# get barcodes for this cell type
    ct_barcodes <- rownames(seurat_obj@meta.data)[
      seurat_obj@meta.data[[broad_col]] == cell_type_name
    ]
    
    if (length(ct_barcodes) < 10) {
      cat("  ✗ Too few cells for radar plot. Skipping.\n")
    } else {
      pdf(file.path(pdir, paste0("P5_RadarPlot_", ct_short, ".pdf")),
          width = 10, height = 8)
      ModuleRadarPlot(
        seurat_obj,
        group.by       = condition_col,
        barcodes       = ct_barcodes,
        axis.label.size = 3.5,
        grid.label.size = 3.5,
        wgcna_name     = wgcna_name
      )
      title(main = paste0(cell_type_name,
                          " – Module Activity Across Conditions\n",
                          "(Radar plot | PCOS and Aging comparison)"),
            cex.main = 1.1, font.main = 2)
      dev.off()
      
      png(file.path(pdir, paste0("P5_RadarPlot_", ct_short, ".png")),
          width = 10, height = 8, units = "in", res = 300)
      ModuleRadarPlot(
        seurat_obj,
        group.by       = condition_col,
        barcodes       = ct_barcodes,
        axis.label.size = 3.5,
        grid.label.size = 3.5,
        wgcna_name     = wgcna_name
      )
      title(main = paste0(cell_type_name,
                          " – Module Activity Across Conditions\n",
                          "(Radar plot | PCOS and Aging comparison)"),
            cex.main = 1.1, font.main = 2)
      dev.off()
      cat("  ✓ Radar plot saved.\n")
    }
  }, error = function(e) {
    cat("  ✗ RadarPlot error:", conditionMessage(e), "\n")
    try(dev.off(), silent = TRUE)
  })
  
# plot p6: modulecorrelogram
# correlation matrix between module eigengenes.
# highly correlated modules may share regulatory programs.
# anti-correlated modules represent opposing biological processes
# (e.g., oxidative phosphorylation vs. inflammatory response).
  cat("[P6] ModuleCorrelogram — module-module correlation...\n")
  tryCatch({
    pdf(file.path(pdir, paste0("P6_ModuleCorrelogram_", ct_short, ".pdf")),
        width = max(8, n_mods * 0.6 + 2),
        height = max(8, n_mods * 0.6 + 2))
    ModuleCorrelogram(
      seurat_obj,
      wgcna_name = wgcna_name
    )
    title(main = paste0(cell_type_name,
                        " – Module Eigengene Correlation Matrix\n",
                        "(Positive = co-regulated | Negative = opposing programs)"),
          cex.main = 1, font.main = 2, line = -1)
    dev.off()
    
    png(file.path(pdir, paste0("P6_ModuleCorrelogram_", ct_short, ".png")),
        width  = max(8, n_mods * 0.6 + 2),
        height = max(8, n_mods * 0.6 + 2),
        units  = "in", res = 300)
    ModuleCorrelogram(
      seurat_obj,
      wgcna_name = wgcna_name
    )
    title(main = paste0(cell_type_name,
                        " – Module Eigengene Correlation Matrix\n",
                        "(Positive = co-regulated | Negative = opposing programs)"),
          cex.main = 1, font.main = 2, line = -1)
    dev.off()
    cat("  ✓ Correlogram saved.\n")
  }, error = function(e) {
    cat("  ✗ Correlogram error:", conditionMessage(e), "\n")
    try(dev.off(), silent = TRUE)
  })
  
# plot p7: dotplot of hmes by condition
# the tutorial's "custom seurat visualization" — adds mes to @meta.data
# then uses dotplot to show average module activity per condition.
# publication-ready summary of which modules are active in which condition.
  cat("[P7] DotPlot — hMEs × condition...\n")
  tryCatch({
    hMEs    <- GetMEs(seurat_obj, harmonized = TRUE, wgcna_name = wgcna_name)
    me_cols <- colnames(hMEs)
    
# add to metadata (overwrite if existing)
    for (col in me_cols) {
      seurat_obj@meta.data[[col]] <- NA_real_
      seurat_obj@meta.data[rownames(hMEs), col] <- hMEs[, col]
    }
    
# subset to this cell type
    ct_cells   <- which(seurat_obj@meta.data[[broad_col]] == cell_type_name)
    seurat_sub <- seurat_obj[, ct_cells]
    
# build dotplot using mes stored in metadata
# we use a custom ggplot because dotplot on meta.data cols is more reliable
    me_long <- seurat_sub@meta.data %>%
      dplyr::select(all_of(c(condition_col, me_cols))) %>%
      pivot_longer(cols = all_of(me_cols),
                   names_to  = "Module",
                   values_to = "ME") %>%
      filter(!is.na(ME)) %>%
      mutate(
        condition = factor(.data[[condition_col]],
                           levels = c("YoungControl","Aged",
                                      "PCOS_Control","PCOS_Case"))
      )
    
    me_summary <- me_long %>%
      group_by(condition, Module) %>%
      summarise(
        avg_ME  = mean(ME, na.rm = TRUE),
        pct_pos = mean(ME > 0, na.rm = TRUE) * 100,
        .groups = "drop"
      )
    
    p_dot <- ggplot(me_summary,
                    aes(x = condition, y = Module,
                        color = avg_ME, size = pct_pos)) +
      geom_point() +
      scale_color_gradient2(
        low      = "#2166AC",
        mid      = "white",
        high     = "#D6604D",
        midpoint = 0,
        name     = "Avg. hME"
      ) +
      scale_size_continuous(range = c(1, 8), name = "% Cells\n(ME > 0)") +
      scale_x_discrete(limits = c("YoungControl","Aged",
                                  "PCOS_Control","PCOS_Case")) +
      labs(
        title    = paste0(cell_type_name,
                          " – Harmonized Module Eigengene Activity by Condition"),
        subtitle = paste0("Dot size = % cells with positive ME | ",
                          "Color = average harmonized ME score\n",
                          "Red = elevated activity | Blue = suppressed activity"),
        x = "Condition",
        y = "Co-expression Module"
      ) +
      theme_cowplot(font_size = 11) +
      theme(
        axis.text.x   = element_text(angle = 35, hjust = 1, size = 10),
        axis.text.y   = element_text(size = 9),
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 12),
        plot.subtitle = element_text(size = 8.5, hjust = 0.5),
        legend.position = "right"
      )
    
    h <- max(5, n_mods * 0.4 + 3)
    ggsave(file.path(pdir, paste0("P7_DotPlot_hME_condition_", ct_short, ".pdf")),
           p_dot, width = 9, height = h, dpi = 300)
    ggsave(file.path(pdir, paste0("P7_DotPlot_hME_condition_", ct_short, ".png")),
           p_dot, width = 9, height = h, dpi = 300)
    cat("  ✓ hME DotPlot saved.\n")
    rm(seurat_sub, me_long, me_summary)
    gc()
  }, error = function(e) {
    cat("  ✗ DotPlot error:", conditionMessage(e), "\n")
  })
  
# table p8: gethubgenes
# the hub gene table is critical for:
# 1. biomarker shortlisting (project objective 4)
# 2. enrichment analysis input
# 3. spatial transcriptomics overlay
  cat("[P8] GetHubGenes — saving top hub genes table...\n")
  tryCatch({
    hub_df <- GetHubGenes(seurat_obj, n_hubs = 20, wgcna_name = wgcna_name)
    hub_df$cell_type <- cell_type_name
    
    write.csv(hub_df,
              file.path(out_dir, paste0("hub_genes_top20_", ct_short, ".csv")),
              row.names = FALSE)
    
# also print top 5 per module to console
    cat("  Top 5 hub genes per module:\n")
    hub_top5 <- hub_df %>%
      group_by(module) %>%
      slice_max(kME, n = 5) %>%
      ungroup()
    print(as.data.frame(hub_top5), row.names = FALSE)
    cat("\n")
    
# hub gene bubble plot — visual version of hub table
    hub_plot_df <- hub_df %>%
      group_by(module) %>%
      slice_max(kME, n = 10) %>%
      ungroup() %>%
      mutate(gene_name = factor(gene_name,
                                levels = rev(unique(gene_name))))
    
    p_hub_bar <- ggplot(hub_plot_df,
                        aes(x = kME, y = gene_name, fill = module)) +
      geom_bar(stat = "identity", width = 0.75) +
      facet_wrap(~module, scales = "free_y", ncol = min(4, n_mods)) +
      scale_fill_manual(
        values = setNames(
          modules$color[match(mods, modules$module)],
          mods
        ),
        guide = "none"
      ) +
      labs(
        title    = paste0(cell_type_name, " – Top Hub Genes by Module (kME)"),
        subtitle = paste0("kME = eigengene-based connectivity | ",
                          "Higher kME = more central hub gene | Top 10 per module"),
        x = "kME (Module Membership)",
        y = "Gene"
      ) +
      theme_cowplot(font_size = 9) +
      theme(
        strip.text    = element_text(face = "bold", size = 9),
        axis.text.y   = element_text(size = 7),
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 12),
        plot.subtitle = element_text(size = 8, hjust = 0.5)
      )
    
    w <- min(20, min(4, n_mods) * 4 + 2)
    h <- ceiling(n_mods / min(4, n_mods)) * 3 + 2
    ggsave(file.path(pdir, paste0("P8_HubGenes_barplot_", ct_short, ".pdf")),
           p_hub_bar, width = w, height = h, dpi = 300, limitsize = FALSE)
    ggsave(file.path(pdir, paste0("P8_HubGenes_barplot_", ct_short, ".png")),
           p_hub_bar, width = w, height = h, dpi = 300, limitsize = FALSE)
    cat("  ✓ Hub genes table and bar plot saved.\n")
  }, error = function(e) {
    cat("  ✗ GetHubGenes error:", conditionMessage(e), "\n")
  })
  
# table p9: getmodules full table
  cat("[P9] GetModules — saving full module assignment table...\n")
  tryCatch({
    mods_full <- GetModules(seurat_obj, wgcna_name = wgcna_name)
    mods_full$cell_type <- cell_type_name
    write.csv(mods_full,
              file.path(out_dir, paste0("module_assignments_", ct_short, ".csv")),
              row.names = FALSE)
    
# module size summary
    mod_sizes <- mods_full %>%
      filter(module != "grey") %>%
      count(module, color) %>%
      arrange(desc(n))
    
    p_sizes <- ggplot(mod_sizes, aes(x = reorder(module, -n), y = n, fill = module)) +
      geom_bar(stat = "identity", width = 0.7) +
      geom_text(aes(label = n), vjust = -0.4, size = 3, fontface = "bold") +
      scale_fill_manual(
        values = setNames(mod_sizes$color, mod_sizes$module),
        guide  = "none"
      ) +
      labs(
        title    = paste0(cell_type_name, " – Co-expression Module Sizes"),
        subtitle = paste0("n = ", n_mods, " modules | Total genes assigned: ",
                          sum(mod_sizes$n)),
        x = "Module",
        y = "Number of Genes"
      ) +
      theme_cowplot(font_size = 11) +
      theme(
        axis.text.x   = element_text(angle = 35, hjust = 1, size = 9),
        plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9)
      )
    
    ggsave(file.path(pdir, paste0("P9_module_sizes_", ct_short, ".pdf")),
           p_sizes, width = max(8, n_mods * 0.7 + 2), height = 5, dpi = 300)
    ggsave(file.path(pdir, paste0("P9_module_sizes_", ct_short, ".png")),
           p_sizes, width = max(8, n_mods * 0.7 + 2), height = 5, dpi = 300)
    cat("  ✓ Module table and size plot saved.\n")
  }, error = function(e) {
    cat("  ✗ GetModules error:", conditionMessage(e), "\n")
  })
  
# plot p10: violin — hme per module × condition
# violin plots show the full distribution of me scores per condition.
# complements the dot plot by showing variance, not just mean.
# wilcoxon test labels indicate statistical significance.
  cat("[P10] Violin plots — hME per module × condition...\n")
  tryCatch({
    hMEs    <- GetMEs(seurat_obj, harmonized = TRUE, wgcna_name = wgcna_name)
    me_cols <- colnames(hMEs)
    
    meta_ct <- seurat_obj@meta.data %>%
      filter(.data[[broad_col]] == cell_type_name) %>%
      dplyr::select(all_of(condition_col))
    
    meta_ct <- cbind(meta_ct,
                     hMEs[rownames(meta_ct), , drop = FALSE])
    
    me_long <- meta_ct %>%
      pivot_longer(cols = all_of(me_cols),
                   names_to  = "Module",
                   values_to = "hME") %>%
      filter(!is.na(hME)) %>%
      mutate(
        Condition = factor(.data[[condition_col]],
                           levels = c("YoungControl","Aged",
                                      "PCOS_Control","PCOS_Case"))
      )
    
# plot top 9 modules if many (to keep readable)
    mods_to_plot <- mods[1:min(9, n_mods)]
    me_filt <- me_long %>% filter(Module %in% mods_to_plot)
    
    p_vln <- ggplot(me_filt,
                    aes(x = Condition, y = hME, fill = Condition)) +
      geom_violin(trim = TRUE, scale = "width", alpha = 0.85) +
      geom_boxplot(width = 0.12, outlier.shape = NA,
                   fill = "white", alpha = 0.7, linewidth = 0.4) +
      geom_hline(yintercept = 0, linetype = "dashed",
                 color = "grey50", linewidth = 0.4) +
      facet_wrap(~Module, scales = "free_y",
                 ncol = min(3, length(mods_to_plot))) +
      scale_fill_manual(values = condition_colors, guide = "none") +
      labs(
        title    = paste0(cell_type_name,
                          " – Module Eigengene Distribution by Condition"),
        subtitle = paste0("Harmonized MEs | Violin + boxplot\n",
                          "Dashed line = ME = 0 (module inactive threshold)"),
        x = "Condition",
        y = "Harmonized ME Score"
      ) +
      theme_cowplot(font_size = 10) +
      theme(
        axis.text.x   = element_text(angle = 35, hjust = 1, size = 8),
        strip.text    = element_text(face = "bold", size = 9),
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 12),
        plot.subtitle = element_text(size = 8.5, hjust = 0.5)
      )
    
    n_mods_plot <- length(mods_to_plot)
    w <- min(3, n_mods_plot) * 4
    h <- ceiling(n_mods_plot / 3) * 4 + 2
    ggsave(file.path(pdir, paste0("P10_violin_hME_condition_", ct_short, ".pdf")),
           p_vln, width = w, height = h, dpi = 300, limitsize = FALSE)
    ggsave(file.path(pdir, paste0("P10_violin_hME_condition_", ct_short, ".png")),
           p_vln, width = w, height = h, dpi = 300, limitsize = FALSE)
    cat("  ✓ ME violin plots saved.\n")
    rm(me_long, me_filt, meta_ct)
    gc()
  }, error = function(e) {
    cat("  ✗ Violin plot error:", conditionMessage(e), "\n")
  })
  
  cat("\n  ✅", cell_type_name, "tutorial plots COMPLETE.\n")
  cat("  Saved to:", pdir, "\n")
  
  return(invisible(seurat_obj))
}

# run for each compartment

# map experiment names to cell type labels
experiments <- list(
  list(wgcna_name = "Granulosa_WGCNA",
       cell_type_name = "Granulosa",
       ct_short = "granulosa"),
  list(wgcna_name = "Stromal_WGCNA",
       cell_type_name = "Stromal",
       ct_short = "stromal"),
  list(wgcna_name = "Immune_WGCNA",
       cell_type_name = "Immune",
       ct_short = "immune")
)

for (exp in experiments) {
  if (exp$wgcna_name %in% names(seurat_obj@misc) ||
      exp$wgcna_name %in% present_wgcna) {
    seurat_obj <- generate_tutorial_plots(
      seurat_obj     = seurat_obj,
      wgcna_name     = exp$wgcna_name,
      cell_type_name = exp$cell_type_name,
      ct_short       = exp$ct_short
    )
    gc()
  } else {
    cat("\n⚠ Skipping", exp$wgcna_name, "— not found in object.\n")
  }
}

# save updated object (with moduleexprscore added)
cat("\nSaving updated object with hub gene scores...\n")
saveRDS(seurat_obj,
        file.path(rds_dir, "seurat_hdwgcna_with_scores.rds"))
cat("✓ Saved:", file.path(rds_dir, "seurat_hdwgcna_with_scores.rds"), "\n")

# done
cat("\n\n", strrep("█", 65), "\n")
cat("  ✅ ALL TUTORIAL PLOTS COMPLETE\n")
cat(strrep("█", 65), "\n\n")
cat("Tutorial plots saved to:", tut_dir, "\n\n")
cat("Files per compartment:\n")
cat("  P1  – PlotDendrogram           (gene dendrogram + module color bar)\n")
cat("  P2  – PlotKMEs                 (kME distributions per module)\n")
cat("  P3  – ModuleFeaturePlot hMEs   (harmonized ME on full UMAP)\n")
cat("  P4  – ModuleFeaturePlot scores (UCell hub gene scores on UMAP)\n")
cat("  P5  – ModuleRadarPlot          (module activity across conditions)\n")
cat("  P6  – ModuleCorrelogram        (module-module correlation matrix)\n")
cat("  P7  – DotPlot hMEs             (module eigengenes × condition)\n")
cat("  P8  – HubGenes bar plot + CSV  (top 20 hub genes per module)\n")
cat("  P9  – Module sizes bar plot    (gene count per module)\n")
cat("  P10 – Violin hME × condition   (full ME distribution)\n\n")

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo_tutorial_plots.txt"))