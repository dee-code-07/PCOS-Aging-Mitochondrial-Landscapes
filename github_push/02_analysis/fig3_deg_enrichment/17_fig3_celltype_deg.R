#!/usr/bin/env Rscript
# task d v3: cell-type-resolved deg — rna assay only
# project: mitochondrial dysfunction in pcos and ovarian aging
# script:  d_celltype_resolved_deg_v3.r
# root cause of v1/v2 failures:
# the annotated objects have sct as default assay. sct findmarkers requires
# prepsctfindmarkers() which fails with >500mb memory error. the previous
# scripts fell through to sct anyway because condition detection failed
# first, so the rna fallback code was never reached.
# this version:
# 1. immediately forces defaultassay = "rna" on load — sct never touched
# 2. joins layers (seurat v5) and runs normalizedata upfront if needed
# 3. adds condition column by pattern-matching sample_id
# 4. runs wilcoxon findmarkers on rna assay — no prepsctfindmarkers needed
# 5. options(future.globals.maxsize = inf) at top to silence future warnings

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(patchwork)
  library(tidyr)
  library(future)
})

# must be set before loading objects
options(future.globals.maxSize = Inf)
plan("sequential")

# configuration

PROJECT_ROOT <- "E:/Documents/mini_project"

PCOS_RDS  <- file.path(PROJECT_ROOT,
                       "scrna/output/07_annotation/pcos/pcos_annotated.rds")
AGING_RDS <- file.path(PROJECT_ROOT,
                       "scrna/output/07_annotation/aging/aging_annotated.rds")
PRIORITY_PATH <- file.path(PROJECT_ROOT,
                           "analysis/09_gene_prioritization/Top_Priority_Genes_Tier1.csv")
OUT_DIR <- file.path(PROJECT_ROOT, "analysis/17_celltype_DEG")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

CELLTYPE_COL <- "final_celltype"

# sample_id patterns → condition groups
# pcos:  "case" → group1,  "ctrl" → group2
# aging: "_a_"  → group1,  "_y_"  → group2
PCOS_CASE_PATTERN   <- "case"
PCOS_CTRL_PATTERN   <- "ctrl"
AGING_OLD_PATTERN   <- "_a_"
AGING_YOUNG_PATTERN <- "_y_"

# deg thresholds
LOG2FC_THRESH <- 0.25
MIN_PCT       <- 0.10
PADJ_THRESH   <- 0.05
MIN_CELLS     <- 10
TOP_N_PLOT    <- 5

# helpers

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, colour = "grey40"),
      axis.text     = element_text(colour = "black"),
      strip.background = element_rect(fill = "grey92", colour = NA),
      strip.text    = element_text(face = "bold")
    )
}

clean_name <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

# load and strip to rna only
# loading only the rna assay saves ram and avoids sct model errors entirely.

load_rna_only <- function(rds_path, dataset_label) {
  
  message("\nLoading ", dataset_label, "...")
  so <- readRDS(rds_path)
  message("  ", ncol(so), " cells | assays: ",
          paste(names(so@assays), collapse = ", "))
  
# force rna as default immediately
  DefaultAssay(so) <- "RNA"
  
# join layers (seurat v5 stores counts/data as separate layers per sample)
  if (inherits(so[["RNA"]], "Assay5")) {
    message("  Joining RNA layers (Seurat v5)...")
    so <- JoinLayers(so, assay = "RNA")
  }
  
# check data layer
  data_sum <- tryCatch(
    sum(GetAssayData(so, assay = "RNA", layer = "data")),
    error = function(e) 0
  )
  message("  RNA data layer sum: ", round(data_sum, 0))
  
  if (data_sum == 0) {
    message("  Data layer is empty — running NormalizeData (LogNormalize)...")
    so <- NormalizeData(so, assay = "RNA",
                        normalization.method = "LogNormalize",
                        scale.factor = 10000, verbose = FALSE)
    message("  NormalizeData complete")
  } else {
    message("  RNA data layer already populated — skipping NormalizeData")
  }
  
  return(so)
}

# add condition column

add_condition <- function(so, dataset) {
  
  sid <- as.character(so@meta.data$sample_id)
  message("  Unique sample_ids: ", paste(sort(unique(sid)), collapse = ", "))
  
  if (dataset == "pcos") {
    cond <- dplyr::case_when(
      grepl(PCOS_CASE_PATTERN, sid, ignore.case = TRUE) ~ "case",
      grepl(PCOS_CTRL_PATTERN, sid, ignore.case = TRUE) ~ "control",
      TRUE ~ NA_character_
    )
    g1 <- "case";    l1 <- "PCOS_case"
    g2 <- "control"; l2 <- "PCOS_ctrl"
  } else {
    cond <- dplyr::case_when(
      grepl(AGING_OLD_PATTERN,   sid, ignore.case = TRUE) ~ "aged",
      grepl(AGING_YOUNG_PATTERN, sid, ignore.case = TRUE) ~ "young",
      TRUE ~ NA_character_
    )
    g1 <- "aged";  l1 <- "Aged"
    g2 <- "young"; l2 <- "Young"
  }
  
  n_na <- sum(is.na(cond))
  if (n_na > 0)
    message("  WARNING: ", n_na, " cells unassigned — check patterns above")
  
  so@meta.data$condition <- cond
  
  tab <- table(cond)
  message("  Condition counts: ",
          paste(names(tab), as.integer(tab), sep = "=", collapse = " | "))
  
  list(so = so, g1 = g1, l1 = l1, g2 = g2, l2 = l2)
}

# run findmarkers per cell type

run_deg <- function(prep, dataset_label) {
  
  so <- prep$so
  g1 <- prep$g1;  l1 <- prep$l1
  g2 <- prep$g2;  l2 <- prep$l2
  
# set condition as active identity
  so@meta.data$condition <- as.character(so@meta.data$condition)
  Idents(so) <- so@meta.data$condition
  
  cell_types <- sort(unique(na.omit(so@meta.data[[CELLTYPE_COL]])))
  message("\n", strrep("-", 60))
  message("  DEG analysis: ", dataset_label, " | ", length(cell_types),
          " cell types | assay: RNA")
  message(strrep("-", 60))
  
  all_deg      <- list()
  summary_rows <- list()
  
  for (ct in cell_types) {
    
    message("\n  ── ", ct, " ──")
    
# subset to this cell type, keeping condition label
    cells_ct <- rownames(so@meta.data)[
      !is.na(so@meta.data[[CELLTYPE_COL]]) &
        so@meta.data[[CELLTYPE_COL]] == ct &
        !is.na(so@meta.data$condition)]
    
    so_ct    <- subset(so, cells = cells_ct)
    Idents(so_ct) <- so_ct@meta.data$condition
    
    n_g1 <- sum(so_ct@meta.data$condition == g1, na.rm = TRUE)
    n_g2 <- sum(so_ct@meta.data$condition == g2, na.rm = TRUE)
    message("  ", l1, "=", n_g1, " | ", l2, "=", n_g2)
    
    if (n_g1 < MIN_CELLS || n_g2 < MIN_CELLS) {
      message("  SKIP (need ≥", MIN_CELLS, " per group)")
      summary_rows[[ct]] <- tibble(
        dataset = dataset_label, cell_type = ct,
        n_g1 = n_g1, n_g2 = n_g2,
        n_up = NA_integer_, n_down = NA_integer_,
        n_sig = NA_integer_, status = "skipped_low_n"
      )
      next
    }
    
# findmarkers — rna assay, log-normalised data layer, wilcoxon
    deg <- tryCatch(
      FindMarkers(
        so_ct,
        ident.1         = g1,
        ident.2         = g2,
        assay           = "RNA",
        layer           = "data",     # explicit: use log-normalised layer
        test.use        = "wilcox",
        logfc.threshold = LOG2FC_THRESH,
        min.pct         = MIN_PCT,
        verbose         = FALSE
      ),
      error = function(e) {
# seurat v4 fallback: uses slot instead of layer
        tryCatch(
          FindMarkers(
            so_ct,
            ident.1         = g1,
            ident.2         = g2,
            assay           = "RNA",
            slot            = "data",
            test.use        = "wilcox",
            logfc.threshold = LOG2FC_THRESH,
            min.pct         = MIN_PCT,
            verbose         = FALSE
          ),
          error = function(e2) {
            message("  ERROR: ", e2$message)
            NULL
          }
        )
      }
    )
    
    if (is.null(deg) || nrow(deg) == 0) {
      message("  No DEGs returned")
      summary_rows[[ct]] <- tibble(
        dataset = dataset_label, cell_type = ct,
        n_g1 = n_g1, n_g2 = n_g2,
        n_up = 0L, n_down = 0L, n_sig = 0L, status = "no_degs"
      )
      next
    }
    
    deg <- deg %>%
      rownames_to_column("gene") %>%
      mutate(
        cell_type   = ct,
        dataset     = dataset_label,
        comparison  = paste0(l1, "_vs_", l2),
        direction   = ifelse(avg_log2FC > 0, "up", "down"),
        significant = p_val_adj < PADJ_THRESH & abs(avg_log2FC) >= LOG2FC_THRESH
      ) %>%
      arrange(p_val_adj, desc(abs(avg_log2FC)))
    
    n_sig  <- sum(deg$significant, na.rm = TRUE)
    n_up   <- sum(deg$significant & deg$direction == "up",   na.rm = TRUE)
    n_down <- sum(deg$significant & deg$direction == "down", na.rm = TRUE)
    message("  Significant: ", n_sig, "  (↑", n_up, "  ↓", n_down, ")")
    
    all_deg[[ct]] <- deg
    summary_rows[[ct]] <- tibble(
      dataset = dataset_label, cell_type = ct,
      n_g1 = n_g1, n_g2 = n_g2,
      n_up = n_up, n_down = n_down,
      n_sig = n_sig, status = "completed"
    )
  }
  
  list(per_celltype = all_deg,
       summary      = bind_rows(summary_rows),
       g1 = g1, l1 = l1, g2 = g2, l2 = l2)
}

# load, prep, run

message("\n", strrep("=", 65))
message("  TASK D v3: CELL-TYPE DEG (RNA ASSAY ONLY)")
message(strrep("=", 65))

so_pcos  <- load_rna_only(PCOS_RDS,  "PCOS")
so_aging <- load_rna_only(AGING_RDS, "Aging")

pcos_prep  <- add_condition(so_pcos,  "pcos")
aging_prep <- add_condition(so_aging, "aging")

# load priority genes
priority_genes <- tryCatch({
  df <- read_csv(PRIORITY_PATH, show_col_types = FALSE)
  gc <- intersect(c("gene","Gene","gene_symbol","Gene_Symbol","symbol"),
                  colnames(df))[1]
  unique(df[[gc]])
}, error = function(e) character(0))
message("\nTier 1 priority genes: ", length(priority_genes))

pcos_res  <- run_deg(pcos_prep,  "PCOS")
aging_res <- run_deg(aging_prep, "Aging")

# save tables

message("\n", strrep("=", 65))
message("  SAVING TABLES")
message(strrep("=", 65))

save_results <- function(res, dataset_label) {
  
# always save the summary
  write_csv(res$summary,
            file.path(OUT_DIR, paste0(dataset_label, "_DEG_summary.csv")))
  
  if (length(res$per_celltype) == 0) {
    message("  No DEG results for ", dataset_label)
    return(invisible(NULL))
  }
  
  out_sub <- file.path(OUT_DIR, tolower(dataset_label))
  dir.create(out_sub, recursive = TRUE, showWarnings = FALSE)
  
  for (ct in names(res$per_celltype)) {
    write_csv(res$per_celltype[[ct]],
              file.path(out_sub,
                        paste0(dataset_label, "_DEG_", clean_name(ct), ".csv")))
  }
  
  combined <- bind_rows(res$per_celltype)
  write_csv(combined,
            file.path(OUT_DIR,
                      paste0(dataset_label, "_DEG_all_celltypes.csv")))
  
  n_done <- sum(res$summary$status == "completed", na.rm = TRUE)
  n_sig  <- sum(res$summary$n_sig,  na.rm = TRUE)
  message("  ✓ ", dataset_label, ": ", n_done,
          " cell types | ", n_sig, " significant DEGs")
}

save_results(pcos_res,  "PCOS")
save_results(aging_res, "Aging")

# priority gene table

priority_deg_df <- NULL

if (length(priority_genes) > 0) {
  rows <- list()
  for (dl in c("PCOS","Aging")) {
    r <- if (dl == "PCOS") pcos_res else aging_res
    for (ct in names(r$per_celltype)) {
      hits <- r$per_celltype[[ct]] %>% filter(gene %in% priority_genes)
      if (nrow(hits) == 0) next
      rows[[length(rows)+1]] <- hits
    }
  }
  if (length(rows) > 0) {
    priority_deg_df <- bind_rows(rows) %>%
      arrange(gene, dataset, cell_type)
    write_csv(priority_deg_df,
              file.path(OUT_DIR, "DEG_priority_gene_celltype.csv"))
    n_sig_pg <- length(unique(
      priority_deg_df$gene[priority_deg_df$significant]))
    message("✓ Priority gene table: ", nrow(priority_deg_df),
            " entries | ", n_sig_pg, " significant")
  }
}

# figure 1: deg summary barplot

all_summ <- bind_rows(pcos_res$summary, aging_res$summary) %>%
  filter(status == "completed")

if (nrow(all_summ) > 0) {
  
  bar_df <- all_summ %>%
    select(dataset, cell_type, n_up, n_down) %>%
    pivot_longer(c(n_up, n_down), names_to = "dir", values_to = "n") %>%
    mutate(
      n         = ifelse(dir == "n_down", -as.numeric(n), as.numeric(n)),
      direction = recode(dir, n_up = "Up", n_down = "Down"),
      cell_type = gsub("_"," ", cell_type)
    )
  
  p_bar <- ggplot(bar_df,
                  aes(x = reorder(cell_type, abs(n), sum),
                      y = n, fill = direction)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_fill_manual(values = c("Up"="#d73027","Down"="#4575b4")) +
    coord_flip() +
    facet_wrap(~dataset, scales = "free_x", ncol = 2) +
    theme_pub() +
    labs(title    = "Significant DEGs per Cell Type",
         subtitle = paste0("FDR<", PADJ_THRESH,
                           "  |  RNA assay (Wilcoxon)"),
         x = NULL, y = "DEGs (up/down)", fill = "Direction")
  
  ggsave(file.path(OUT_DIR, "DEG_summary_barplot.png"),
         p_bar, width = 12, height = 6, dpi = 300)
  ggsave(file.path(OUT_DIR, "DEG_summary_barplot.tiff"),
         p_bar, width = 12, height = 6, dpi = 600, compression = "lzw")
  message("✓ DEG summary barplot saved")
}

# figure 2: dot plot

plot_dotplot <- function(res, so, dataset_label) {
  
  if (length(res$per_celltype) == 0) return(invisible(NULL))
  
  top_df <- bind_rows(res$per_celltype) %>%
    filter(significant) %>%
    group_by(cell_type) %>%
    slice_max(order_by = abs(avg_log2FC), n = TOP_N_PLOT) %>%
    ungroup()
  
  if (nrow(top_df) == 0) {
    message("  No significant DEGs for dotplot (", dataset_label, ")")
    return(invisible(NULL))
  }
  
  genes_plot <- unique(top_df$gene)
  message("  Dotplot: ", length(genes_plot), " genes (", dataset_label, ")")
  
  Idents(so) <- so@meta.data[[CELLTYPE_COL]]
  DefaultAssay(so) <- "RNA"
  
  p <- tryCatch(
    DotPlot(so, features = genes_plot, assay = "RNA",
            cols = c("lightgrey","#d73027"), dot.scale = 5) +
      coord_flip() +
      theme_pub(base_size = 9) +
      theme(axis.text.x = element_text(angle=45, hjust=1, size=7.5),
            axis.text.y = element_text(size=7.5, face="italic")) +
      labs(title    = paste0(dataset_label, " — Top DEGs per Cell Type"),
           subtitle = paste0("Top ", TOP_N_PLOT,
                             " DEGs by |log₂FC| | RNA assay"),
           x = "Gene", y = "Cell Type"),
    error = function(e) { message("  DotPlot error: ", e$message); NULL }
  )
  
  if (!is.null(p)) {
    h <- max(5, 0.25 * length(genes_plot) + 2)
    w <- max(8, 0.6  * length(unique(so@meta.data[[CELLTYPE_COL]])) + 3)
    ggsave(file.path(OUT_DIR, paste0(dataset_label, "_DEG_dotplot.png")),
           p, width = w, height = h, dpi = 300)
    ggsave(file.path(OUT_DIR, paste0(dataset_label, "_DEG_dotplot.tiff")),
           p, width = w, height = h, dpi = 600, compression = "lzw")
    message("✓ Dotplot saved: ", dataset_label)
  }
}

plot_dotplot(pcos_res,  pcos_prep$so,  "PCOS")
plot_dotplot(aging_res, aging_prep$so, "Aging")

# figure 3: priority gene fc heatmap

if (!is.null(priority_deg_df)) {
  
  sig_df <- priority_deg_df %>%
    filter(significant) %>%
    mutate(col_label = paste0(dataset, "\n", gsub("_"," ", cell_type)))
  
  if (nrow(sig_df) > 0) {
    fc_wide <- sig_df %>%
      select(gene, col_label, avg_log2FC) %>%
      pivot_wider(names_from  = col_label, values_from = avg_log2FC,
                  values_fn   = mean, values_fill = 0)
    
    fc_mat <- as.matrix(fc_wide[,-1])
    rownames(fc_mat) <- fc_wide$gene
    fc_mat <- fc_mat[order(rowMeans(fc_mat), decreasing=TRUE),,drop=FALSE]
    
    fc_long <- as.data.frame(fc_mat) %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to="group", values_to="log2FC") %>%
      mutate(gene  = factor(gene,  levels=rev(rownames(fc_mat))),
             group = factor(group, levels=colnames(fc_mat)))
    
    max_fc <- max(abs(fc_long$log2FC), na.rm=TRUE)
    
    p_heat <- ggplot(fc_long, aes(x=group, y=gene, fill=log2FC)) +
      geom_tile(colour="white", linewidth=0.3) +
      scale_fill_gradient2(low="#4575b4", mid="white", high="#d73027",
                           midpoint=0, limits=c(-max_fc, max_fc),
                           name="log₂FC") +
      theme_pub(base_size=9) +
      theme(axis.text.x = element_text(angle=45, hjust=1, size=8),
            axis.text.y = element_text(face="italic", size=8),
            panel.grid  = element_blank(),
            axis.ticks  = element_blank()) +
      labs(title    = "Tier 1 Priority Genes — log₂FC Across Cell Types",
           subtitle = paste0("Significant DEGs only (FDR<", PADJ_THRESH,")"),
           x=NULL, y=NULL)
    
    h <- max(4, 0.35*nrow(fc_mat)+1.5)
    w <- max(7, 0.55*ncol(fc_mat)+2)
    ggsave(file.path(OUT_DIR,"Priority_DEG_heatmap.png"),
           p_heat, width=w, height=h, dpi=300)
    ggsave(file.path(OUT_DIR,"Priority_DEG_heatmap.tiff"),
           p_heat, width=w, height=h, dpi=600, compression="lzw")
    message("✓ Priority gene heatmap saved")
  }
}

# final summary

message("\n", strrep("=", 65))
message("  TASK D COMPLETE")
message(strrep("=", 65))

for (dl in c("PCOS","Aging")) {
  r <- if (dl=="PCOS") pcos_res else aging_res
  message("\n  ", dl, " DEG Summary:")
  print(as.data.frame(r$summary))
}

message("\nFiles in output directory:")
for (f in list.files(OUT_DIR, recursive=TRUE, full.names=FALSE))
  message("  ", f)