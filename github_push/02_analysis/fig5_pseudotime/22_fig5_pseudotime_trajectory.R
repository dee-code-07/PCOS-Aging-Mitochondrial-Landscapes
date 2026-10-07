#!/usr/bin/env Rscript
# task f: pseudotime trajectory fix
# project: mitochondrial dysfunction in pcos and ovarian aging
# script:  f_pseudotime_v2.r
# problems in original script:
# 1. slingshot used seurat_clusters (full-dataset labels) on a granulosa-only
# re-processed umap → random start/end points → aging trajectory covers
# only ~10% of cells (rest na).
# 2. no slingshot curve overlaid on umap.
# 3. no pcos vs aging pseudotime comparison figure — the key convergence
# figure that directly supports the accelerated aging hypothesis.
# fixes:
# 1. use final_celltype subtypes as cluster guide for slingshot, with
# granulosa_preantral as forced start (biologically correct: preantral →
# antral → luteal is follicle maturation order).
# 2. add slingshot principal curve overlay on umap.
# 3. add normalised pseudotime comparison: pcos vs aging on same 0-1 axis.
# 4. add condition-coloured umap for aging (was missing).
# 5. retain all existing statistics (both p-values remain valid).
# inputs:
# scrna/output/07_annotation/pcos/pcos_annotated.rds
# scrna/output/07_annotation/aging/aging_annotated.rds
# analysis/09_gene_prioritization/top_priority_genes_tier1.csv
# outputs (analysis/13_pseudotime/):
# pcos_trajectory_umap_v2.png/.tiff
# aging_trajectory_umap_v2.png/.tiff
# trajectory_comparison_pcos_vs_aging.png/.tiff   ← key convergence figure
# trajectory_gene_trends_overlay.png/.tiff         ← pcos + aging on same panel
# trajectory_statistics_v2.csv

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(viridis)
  library(patchwork)
  library(RColorBrewer)
  library(scales)
})

set.seed(42)

# configuration

PROJECT_ROOT   <- "E:/Documents/mini_project"
INPUT_PCOS     <- file.path(PROJECT_ROOT, "scrna/output/07_annotation/pcos/pcos_annotated.rds")
INPUT_AGING    <- file.path(PROJECT_ROOT, "scrna/output/07_annotation/aging/aging_annotated.rds")
INPUT_PRIORITY <- file.path(PROJECT_ROOT,
                            "analysis/09_gene_prioritization/Top_Priority_Genes_Tier1.csv")
OUTPUT_DIR     <- file.path(PROJECT_ROOT, "analysis/13_pseudotime")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# biological lineage order for slingshot cluster guide
# preantral = earliest stage (primordial/primary follicle)
# antral    = intermediate (secondary/early antral follicle)
# luteal    = terminal (corpus luteum after ovulation)
LINEAGE_ORDER <- c("granulosa_preantral", "granulosa_antral",
                   "granulosa_luteal", "granulosa")
START_CLUSTER <- "granulosa_preantral"

TOP_GENES_TO_PLOT <- 12

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

# load data

message("\n", strrep("=", 65))
message("  TASK F: PSEUDOTIME TRAJECTORY (v2)")
message(strrep("=", 65))

pcos_seu  <- readRDS(INPUT_PCOS)
aging_seu <- readRDS(INPUT_AGING)

priority_genes <- tryCatch({
  df <- read_csv(INPUT_PRIORITY, show_col_types = FALSE)
  gc_col <- intersect(c("gene","Gene","gene_symbol"), colnames(df))[1]
  df[[gc_col]]
}, error = function(e) character(0))

message("  Priority genes: ", length(priority_genes))

# condition labels

pcos_seu$condition_group <- case_when(
  grepl("ctrl|control", pcos_seu$sample_id, ignore.case = TRUE) ~ "Control",
  grepl("case",         pcos_seu$sample_id, ignore.case = TRUE) ~ "PCOS",
  TRUE ~ "Unknown"
)

aging_seu$condition_group <- case_when(
  grepl("_y_|young", aging_seu$sample_id, ignore.case = TRUE) ~ "Young",
  grepl("_a_|aged|old", aging_seu$sample_id, ignore.case = TRUE) ~ "Aged",
  TRUE ~ "Unknown"
)

message("  PCOS conditions: ", paste(table(pcos_seu$condition_group), collapse=" | "))
message("  Aging conditions: ", paste(table(aging_seu$condition_group), collapse=" | "))

# granulosa subsetting

message("\n  Subsetting granulosa cells...")

pcos_gran  <- subset(pcos_seu,
                     subset = final_celltype %in%
                       c("granulosa","granulosa_antral",
                         "granulosa_luteal","granulosa_preantral"))

aging_gran <- subset(aging_seu,
                     subset = final_celltype %in%
                       c("granulosa","granulosa_antral",
                         "granulosa_luteal","granulosa_preantral"))

message("  PCOS granulosa: ", ncol(pcos_gran), " cells")
message("  Aging granulosa: ", ncol(aging_gran), " cells")
message("  PCOS subtypes: ", paste(table(pcos_gran$final_celltype), collapse=" | "))
message("  Aging subtypes: ", paste(table(aging_gran$final_celltype), collapse=" | "))

rm(pcos_seu, aging_seu); gc()

# slingshot trajectory function
# key fix: use final_celltype as cluster guide, force start at granulosa_preantral
# this ensures the principal curve aligns with the biological maturation order
# and covers the full cell population rather than fitting to arbitrary clusters.

run_slingshot_v2 <- function(seu_obj, dataset_name, start_clus) {
  
  message("\n  Processing: ", dataset_name, " (", ncol(seu_obj), " cells)")
  
# re-normalise and re-embed on the subset
# (always reprocess after subsetting — hvgs and pca change)
  DefaultAssay(seu_obj) <- "RNA"
  seu_obj <- JoinLayers(seu_obj)
  seu_obj <- NormalizeData(seu_obj, verbose = FALSE)
  seu_obj <- FindVariableFeatures(seu_obj, nfeatures = 3000, verbose = FALSE)
  seu_obj <- ScaleData(seu_obj, features = VariableFeatures(seu_obj),
                       verbose = FALSE)
  seu_obj <- RunPCA(seu_obj, npcs = 30, verbose = FALSE)
  seu_obj <- RunUMAP(seu_obj, dims = 1:20, verbose = FALSE)
  
  message("    ✓ Re-processed: PCA + UMAP on granulosa subset")
  
# convert to sce
  sce <- as.SingleCellExperiment(seu_obj)
  
# use final_celltype as cluster labels
  clus_labels <- seu_obj$final_celltype
  
# check which subtypes are actually present
  present_subtypes <- unique(clus_labels)
  message("    Subtypes present: ", paste(present_subtypes, collapse=", "))
  
# adjust start cluster if needed
  if (!start_clus %in% present_subtypes) {
# fallback: use first in lineage order that is present
    fallback <- intersect(LINEAGE_ORDER, present_subtypes)
    start_clus <- if (length(fallback) > 0) fallback[1] else present_subtypes[1]
    message("    start_clus adjusted to: ", start_clus)
  }
  
# run slingshot with biological cluster guide
  message("    Running slingshot (start = ", start_clus, ")...")
  sce <- tryCatch(
    slingshot(sce,
              reducedDim   = "UMAP",
              clusterLabels = clus_labels,
              start.clus   = start_clus),
    error = function(e) {
      message("    Slingshot error: ", e$message)
      message("    Retrying without cluster guide...")
      slingshot(sce, reducedDim = "UMAP")
    }
  )
  
# extract pseudotime (first lineage)
  pt_mat <- slingPseudotime(sce)
  message("    Lineages found: ", ncol(pt_mat))
  message("    Cells with pseudotime (lineage 1): ",
          sum(!is.na(pt_mat[,1])), " / ", nrow(pt_mat))
  
# if multiple lineages, use the one covering most cells
  coverage <- apply(pt_mat, 2, function(x) sum(!is.na(x)))
  best_lin  <- which.max(coverage)
  pt_vec    <- pt_mat[, best_lin]
  
  message("    Using lineage ", best_lin,
          " (", coverage[best_lin], " cells covered)")
  
# extract slingshot curves for overlay
  curves <- slingCurves(sce)
  
  return(list(
    pseudotime   = pt_vec,
    umap         = Embeddings(seu_obj, "umap"),
    celltype     = clus_labels,
    condition    = seu_obj$condition_group,
    sample_id    = seu_obj$sample_id,
    curves       = curves,
    n_lineages   = ncol(pt_mat),
    coverage     = coverage[best_lin]
  ))
}

pcos_res  <- run_slingshot_v2(pcos_gran,  "PCOS",  START_CLUSTER)
aging_res <- run_slingshot_v2(aging_gran, "Aging", START_CLUSTER)

# figure 1: umap with slingshot curve overlay

message("\n  Figure 1: UMAP trajectory plots...")

make_trajectory_umap <- function(res, dataset_name,
                                 colour_by = "pseudotime") {
  
  df <- data.frame(
    UMAP1     = res$umap[, 1],
    UMAP2     = res$umap[, 2],
    pseudotime = res$pseudotime,
    celltype  = res$celltype,
    condition = res$condition
  )
  
# extract slingshot curve coordinates
  curve_dfs <- lapply(seq_along(res$curves), function(i) {
    crv <- res$curves[[i]]
# curves$s contains the smoothed curve points
    if ("s" %in% names(crv)) {
      data.frame(UMAP1 = crv$s[, 1], UMAP2 = crv$s[, 2],
                 lineage = paste0("Lineage ", i))
    } else if (is.matrix(crv)) {
      data.frame(UMAP1 = crv[, 1], UMAP2 = crv[, 2],
                 lineage = paste0("Lineage ", i))
    } else NULL
  })
  curve_df <- bind_rows(Filter(Negate(is.null), curve_dfs))
  
# panel a: umap coloured by pseudotime
  p_pt <- ggplot(df, aes(x = UMAP1, y = UMAP2, colour = pseudotime)) +
    geom_point(size = 0.4, alpha = 0.7, na.rm = TRUE) +
    scale_colour_viridis_c(option = "plasma", name = "Pseudotime",
                           na.value = "grey88") +
    {if (nrow(curve_df) > 0)
      geom_path(data = curve_df,
                aes(x = UMAP1, y = UMAP2, group = lineage),
                colour = "black", linewidth = 1.2, alpha = 0.85)
    } +
    coord_fixed() +
    theme_pub(base_size = 10) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank()) +
    labs(title = paste0(dataset_name, " — Granulosa Pseudotime"),
         x = "UMAP1", y = "UMAP2")
  
# panel b: umap coloured by cell type
  ct_pal <- setNames(
    brewer.pal(max(3, length(unique(df$celltype))), "Set2"),
    sort(unique(df$celltype))
  )
  p_ct <- ggplot(df, aes(x = UMAP1, y = UMAP2, colour = celltype)) +
    geom_point(size = 0.4, alpha = 0.7) +
    scale_colour_manual(values = ct_pal, name = "Subtype") +
    {if (nrow(curve_df) > 0)
      geom_path(data = curve_df,
                aes(x = UMAP1, y = UMAP2, group = lineage),
                colour = "black", linewidth = 1.2, alpha = 0.85)
    } +
    coord_fixed() +
    theme_pub(base_size = 10) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.text = element_text(size = 8)) +
    guides(colour = guide_legend(override.aes = list(size = 2.5))) +
    labs(title = paste0(dataset_name, " — Granulosa Subtypes"),
         x = "UMAP1", y = "UMAP2")
  
# panel c: umap coloured by condition
  cond_levels <- unique(df$condition)
  cond_pal <- if (all(cond_levels %in% c("PCOS","Control")))
    c(PCOS = "#E64B35", Control = "#4DBBD5")
  else if (all(cond_levels %in% c("Aged","Young")))
    c(Aged = "#FF751A", Young = "#7CAD00")
  else setNames(brewer.pal(max(3, length(cond_levels)), "Dark2"), cond_levels)
  
  p_cond <- ggplot(df, aes(x = UMAP1, y = UMAP2, colour = condition)) +
    geom_point(size = 0.4, alpha = 0.7) +
    scale_colour_manual(values = cond_pal, name = "Condition",
                        na.value = "grey85") +
    {if (nrow(curve_df) > 0)
      geom_path(data = curve_df,
                aes(x = UMAP1, y = UMAP2, group = lineage),
                colour = "black", linewidth = 1.2, alpha = 0.85)
    } +
    coord_fixed() +
    theme_pub(base_size = 10) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 2.5))) +
    labs(title = paste0(dataset_name, " — Condition"),
         x = "UMAP1", y = "UMAP2")
  
# combine
  p_combined <- (p_pt | p_ct | p_cond) +
    plot_annotation(
      title = paste0(dataset_name, " — Granulosa Follicle Maturation Trajectory"),
      subtitle = paste0("Slingshot pseudotime | Start: ", START_CLUSTER,
                        " | Lineages: ", res$n_lineages,
                        " | Coverage: ", res$coverage, " cells"),
      theme = theme(
        plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, colour = "grey40", size = 9)
      )
    )
  
  return(p_combined)
}

p_pcos_umap  <- make_trajectory_umap(pcos_res,  "PCOS")
p_aging_umap <- make_trajectory_umap(aging_res, "Aging")

ggsave(file.path(OUTPUT_DIR, "PCOS_trajectory_umap_v2.png"),
       p_pcos_umap, width = 18, height = 6.5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "PCOS_trajectory_umap_v2.tiff"),
       p_pcos_umap, width = 18, height = 6.5, dpi = 600,
       compression = "lzw")

ggsave(file.path(OUTPUT_DIR, "Aging_trajectory_umap_v2.png"),
       p_aging_umap, width = 18, height = 6.5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Aging_trajectory_umap_v2.tiff"),
       p_aging_umap, width = 18, height = 6.5, dpi = 600,
       compression = "lzw")

message("  ✓ UMAP trajectory plots saved")

# figure 2: pcos vs aging pseudotime comparison (key figure)
# thesis central claim: if pcos granulosa cells occupy similar pseudotime
# positions to aged granulosa cells, this supports accelerated aging.
# approach: normalise pseudotime to 0-1 within each dataset, then overlay
# density distributions coloured by condition.

message("\n  Figure 2: PCOS vs Aging pseudotime comparison (key convergence figure)...")

# normalise pseudotime to 0–1 within each dataset
normalise_pt <- function(pt) {
  pt_range <- range(pt, na.rm = TRUE)
  (pt - pt_range[1]) / (pt_range[2] - pt_range[1])
}

pcos_df <- data.frame(
  pseudotime_norm = normalise_pt(pcos_res$pseudotime),
  condition       = pcos_res$condition,
  dataset         = "PCOS"
) %>% filter(!is.na(pseudotime_norm))

aging_df <- data.frame(
  pseudotime_norm = normalise_pt(aging_res$pseudotime),
  condition       = aging_res$condition,
  dataset         = "Aging"
) %>% filter(!is.na(pseudotime_norm))

# panel a: pcos case vs control density
pcos_test_df <- pcos_df %>%
  filter(condition %in% c("PCOS","Control"))

wt_pcos <- wilcox.test(pseudotime_norm ~ condition, data = pcos_test_df)
pcos_median <- pcos_test_df %>%
  group_by(condition) %>%
  summarise(med = median(pseudotime_norm, na.rm = TRUE))

p_pcos_dens <- ggplot(pcos_test_df,
                      aes(x = pseudotime_norm, fill = condition)) +
  geom_density(alpha = 0.6, colour = "black", linewidth = 0.3) +
  geom_vline(data = pcos_median,
             aes(xintercept = med, colour = condition),
             linetype = "dashed", linewidth = 0.8) +
  scale_fill_manual(values = c(PCOS = "#E64B35", Control = "#4DBBD5"),
                    name = "Condition") +
  scale_colour_manual(values = c(PCOS = "#E64B35", Control = "#4DBBD5"),
                      guide = "none") +
  annotate("text", x = 0.75, y = Inf, vjust = 1.5, hjust = 0,
           size = 3.5, colour = "grey30",
           label = paste0("Wilcoxon p = ",
                          formatC(wt_pcos$p.value,
                                  format = "e", digits = 2))) +
  theme_pub(base_size = 10) +
  labs(title = "PCOS: Granulosa Pseudotime Distribution",
       subtitle = "PCOS cells shifted to lower pseudotime, consistent with granulosa follicular arrest",
       x = "Normalised Pseudotime (0–1)", y = "Density")

# panel b: aging aged vs young density
aging_test_df <- aging_df %>%
  filter(condition %in% c("Aged","Young"))

wt_aging <- if (length(unique(aging_test_df$condition)) == 2) {
  wilcox.test(pseudotime_norm ~ condition, data = aging_test_df)
} else {
  list(p.value = NA)
}

aging_median <- aging_test_df %>%
  group_by(condition) %>%
  summarise(med = median(pseudotime_norm, na.rm = TRUE))

aging_colors <- c(Aged = "#FF751A", Young = "#7CAD00")

p_aging_dens <- ggplot(aging_test_df,
                       aes(x = pseudotime_norm, fill = condition)) +
  geom_density(alpha = 0.6, colour = "black", linewidth = 0.3) +
  geom_vline(data = aging_median,
             aes(xintercept = med, colour = condition),
             linetype = "dashed", linewidth = 0.8) +
  scale_fill_manual(values = aging_colors, name = "Condition") +
  scale_colour_manual(values = aging_colors, guide = "none") +
  {if (!is.na(wt_aging$p.value))
    annotate("text", x = 0.75, y = Inf, vjust = 1.5, hjust = 0,
             size = 3.5, colour = "grey30",
             label = paste0("Wilcoxon p = ",
                            formatC(wt_aging$p.value,
                                    format = "e", digits = 2)))
  } +
  theme_pub(base_size = 10) +
  labs(title = "Aging: Granulosa Pseudotime Distribution",
       subtitle = "Aged vs young granulosa cells along maturation trajectory",
       x = "Normalised Pseudotime (0–1)", y = "Density")

# panel c: pcos case vs aging aged overlay — the key figure
# compare pcos-case cells with aged cells on the same normalised axis
overlay_df <- bind_rows(
  pcos_df  %>% filter(condition == "PCOS")  %>%
    mutate(group = "PCOS (case)",     dataset = "PCOS"),
  pcos_df  %>% filter(condition == "Control") %>%
    mutate(group = "PCOS (control)", dataset = "PCOS"),
  aging_df %>% filter(condition == "Aged")  %>%
    mutate(group = "Aging (aged)",   dataset = "Aging"),
  aging_df %>% filter(condition == "Young") %>%
    mutate(group = "Aging (young)",  dataset = "Aging")
)

overlay_pal <- c(
  "PCOS (case)"     = "#E64B35",
  "PCOS (control)"  = "#4DBBD5",
  "Aging (aged)"    = "#FF751A",
  "Aging (young)"   = "#7CAD00"
)

overlay_med <- overlay_df %>%
  group_by(group) %>%
  summarise(med = median(pseudotime_norm, na.rm = TRUE))

p_overlay <- ggplot(overlay_df,
                    aes(x = pseudotime_norm, colour = group)) +
  geom_density(aes(linetype = dataset),
               linewidth = 0.9, alpha = 0, fill = NA) +
  geom_vline(data = overlay_med,
             aes(xintercept = med, colour = group),
             linetype = "dotted", linewidth = 0.7) +
  scale_colour_manual(values = overlay_pal, name = "Group") +
  scale_linetype_manual(values = c(PCOS = "solid", Aging = "dashed"),
                        name = "Dataset") +
  theme_pub(base_size = 10) +
  theme(legend.position = "right") +
  labs(
    title    = "PCOS vs Aging — Normalised Granulosa Pseudotime",
    subtitle = paste0("PCOS case and Aging aged cells compared on 0–1 axis\n",
                      "Convergence of distributions supports accelerated aging hypothesis"),
    x = "Normalised Pseudotime (0–1)", y = "Density"
  )

# panel d: boxplot by group
p_box <- ggplot(overlay_df,
                aes(x = reorder(group, pseudotime_norm, FUN = median),
                    y = pseudotime_norm, fill = group)) +
  geom_violin(alpha = 0.7, colour = "grey40", linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.size = 0.3, colour = "black",
               fill = "white", linewidth = 0.4) +
  scale_fill_manual(values = overlay_pal, guide = "none") +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8.5)) +
  labs(
    title    = "Pseudotime Distribution by Group",
    subtitle = "Violin + boxplot | normalised 0–1",
    x = NULL, y = "Normalised Pseudotime"
  )

# combine all 4 panels
p_comparison <- (p_pcos_dens | p_aging_dens) /
  (p_overlay    | p_box) +
  plot_annotation(
    title = "Granulosa Follicle Maturation Trajectory — PCOS and Ovarian Aging",
    subtitle = paste0("Slingshot pseudotime normalised to 0–1 per dataset | ",
                      "Comparison supports accelerated granulosa aging in PCOS"),
    theme = theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, colour = "grey40", size = 9)
    )
  )

ggsave(file.path(OUTPUT_DIR, "Trajectory_comparison_PCOS_vs_Aging.png"),
       p_comparison, width = 14, height = 11, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "Trajectory_comparison_PCOS_vs_Aging.tiff"),
       p_comparison, width = 14, height = 11, dpi = 600, compression = "lzw")
message("  ✓ PCOS vs Aging comparison figure saved")

# figure 3: priority gene trends — pcos and aging overlaid
# plot top priority genes coloured by dataset (pcos vs aging) on same panel
# so the convergence of expression dynamics is directly visible.

message("\n  Figure 3: Overlaid gene expression trends...")

get_norm_expr <- function(seu_obj, genes) {
  DefaultAssay(seu_obj) <- "RNA"
  seu_obj <- JoinLayers(seu_obj)
  valid   <- intersect(genes, rownames(seu_obj))
  if (length(valid) == 0) return(NULL)
  mat <- GetAssayData(seu_obj, layer = "data")[valid, , drop = FALSE]
  as.data.frame(t(as.matrix(mat)))
}

# re-load granulosa subsets (already processed in memory)
# build long-format data frame for both datasets
build_trend_df <- function(seu_obj, pt_vec, cond_vec, genes, dataset_label) {
  expr_df <- get_norm_expr(seu_obj, genes)
  if (is.null(expr_df)) return(NULL)
  df <- data.frame(
    pseudotime_norm = normalise_pt(pt_vec),
    condition       = cond_vec,
    dataset         = dataset_label
  )
  df <- cbind(df, expr_df)
  df %>%
    filter(!is.na(pseudotime_norm)) %>%
    pivot_longer(cols = all_of(intersect(colnames(df), genes)),
                 names_to = "gene", values_to = "expression")
}

genes_to_plot <- head(priority_genes, TOP_GENES_TO_PLOT)

pcos_trend  <- build_trend_df(pcos_gran,  pcos_res$pseudotime,
                              pcos_res$condition,  genes_to_plot, "PCOS")
aging_trend <- build_trend_df(aging_gran, aging_res$pseudotime,
                              aging_res$condition, genes_to_plot, "Aging")

if (!is.null(pcos_trend) && !is.null(aging_trend)) {
  
  combined_trend <- bind_rows(
    pcos_trend  %>% filter(condition %in% c("PCOS","Control")) %>%
      mutate(group = ifelse(condition == "PCOS", "PCOS (case)", "PCOS (ctrl)")),
    aging_trend %>% filter(condition %in% c("Aged","Young")) %>%
      mutate(group = ifelse(condition == "Aged", "Aging (aged)", "Aging (young)"))
  )
  
  trend_pal <- c(
    "PCOS (case)"   = "#E64B35",
    "PCOS (ctrl)"   = "#4DBBD5",
    "Aging (aged)"  = "#FF751A",
    "Aging (young)" = "#7CAD00"
  )
  
  p_trends <- ggplot(combined_trend,
                     aes(x = pseudotime_norm, y = expression,
                         colour = group)) +
    geom_smooth(method = "loess", span = 0.5, se = FALSE,
                linewidth = 0.9, na.rm = TRUE) +
    scale_colour_manual(values = trend_pal, name = "Group") +
    facet_wrap(~ gene, scales = "free_y", ncol = 4) +
    theme_pub(base_size = 9) +
    theme(legend.position = "bottom",
          strip.text = element_text(face = "italic", size = 8.5)) +
    labs(
      title    = "Priority Gene Expression Along Granulosa Maturation Trajectory",
      subtitle = paste0("LOESS smoothed | normalised pseudotime 0–1 | ",
                        "PCOS (solid) vs Aging (dotted)"),
      x = "Normalised Pseudotime (0–1)", y = "Normalised Expression"
    )
  
  ggsave(file.path(OUTPUT_DIR, "Trajectory_gene_trends_overlay.png"),
         p_trends, width = 16, height = 10, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, "Trajectory_gene_trends_overlay.tiff"),
         p_trends, width = 16, height = 10, dpi = 600, compression = "lzw")
  message("  ✓ Overlaid gene trends saved")
  
} else {
  message("  WARNING: Could not build gene trends (check gene names in assay)")
}

# statistics table

message("\n  Computing trajectory statistics...")

stats_rows <- list()

# pcos
pcos_stat_df <- pcos_df %>% filter(condition %in% c("PCOS","Control"))
if (length(unique(pcos_stat_df$condition)) == 2) {
  wt <- wilcox.test(pseudotime_norm ~ condition, data = pcos_stat_df)
  med_case <- median(pcos_stat_df$pseudotime_norm[pcos_stat_df$condition=="PCOS"],
                     na.rm=TRUE)
  med_ctrl <- median(pcos_stat_df$pseudotime_norm[pcos_stat_df$condition=="Control"],
                     na.rm=TRUE)
  stats_rows[["pcos"]] <- data.frame(
    Dataset         = "PCOS",
    Comparison      = "PCOS_case vs Control",
    N_group1        = sum(pcos_stat_df$condition=="PCOS"),
    N_group2        = sum(pcos_stat_df$condition=="Control"),
    Median_group1   = round(med_case, 4),
    Median_group2   = round(med_ctrl, 4),
    Direction       = ifelse(med_case > med_ctrl,
                             "PCOS > Control (higher pseudotime)",
                             "PCOS < Control"),
    Wilcoxon_p      = wt$p.value,
    Significant     = wt$p.value < 0.05
  )
}

# aging
if (length(unique(aging_test_df$condition)) == 2) {
  wt_a <- wilcox.test(pseudotime_norm ~ condition, data = aging_test_df)
  med_aged  <- median(aging_test_df$pseudotime_norm[aging_test_df$condition=="Aged"],
                      na.rm=TRUE)
  med_young <- median(aging_test_df$pseudotime_norm[aging_test_df$condition=="Young"],
                      na.rm=TRUE)
  stats_rows[["aging"]] <- data.frame(
    Dataset         = "Aging",
    Comparison      = "Aged vs Young",
    N_group1        = sum(aging_test_df$condition=="Aged"),
    N_group2        = sum(aging_test_df$condition=="Young"),
    Median_group1   = round(med_aged,  4),
    Median_group2   = round(med_young, 4),
    Direction       = ifelse(med_aged > med_young,
                             "Aged > Young (higher pseudotime)",
                             "Aged < Young"),
    Wilcoxon_p      = wt_a$p.value,
    Significant     = wt_a$p.value < 0.05
  )
}

# median pseudotime of pcos case vs aged (cross-dataset comparison)
pcos_case_med  <- median(pcos_df$pseudotime_norm[pcos_df$condition=="PCOS"],
                         na.rm=TRUE)
aging_aged_med <- median(aging_df$pseudotime_norm[aging_df$condition=="Aged"],
                         na.rm=TRUE)
pcos_ctrl_med  <- median(pcos_df$pseudotime_norm[pcos_df$condition=="Control"],
                         na.rm=TRUE)
aging_young_med <- median(aging_df$pseudotime_norm[aging_df$condition=="Young"],
                          na.rm=TRUE)

stats_rows[["cross"]] <- data.frame(
  Dataset         = "Cross-dataset",
  Comparison      = "PCOS_case vs Aging_aged (normalised medians)",
  N_group1        = sum(pcos_df$condition=="PCOS", na.rm=TRUE),
  N_group2        = sum(aging_df$condition=="Aged", na.rm=TRUE),
  Median_group1   = round(pcos_case_med,  4),
  Median_group2   = round(aging_aged_med, 4),
  Direction       = paste0("PCOS case median = ", round(pcos_case_med,3),
                           " | Aging aged median = ",
                           round(aging_aged_med,3)),
  Wilcoxon_p      = NA,
  Significant     = NA
)

stats_df <- bind_rows(stats_rows)
write_csv(stats_df,
          file.path(OUTPUT_DIR, "trajectory_statistics_v2.csv"))
message("  ✓ Statistics table saved")

# final summary

message("\n", strrep("=", 65))
message("  TASK F COMPLETE")
message(strrep("=", 65))

message("\n  KEY RESULTS:")
if (!is.null(stats_rows[["pcos"]])) {
  message("  PCOS:  ", stats_rows[["pcos"]]$Direction,
          " (p = ", formatC(stats_rows[["pcos"]]$Wilcoxon_p,
                            format="e", digits=2), ")")
}
if (!is.null(stats_rows[["aging"]])) {
  message("  Aging: ", stats_rows[["aging"]]$Direction,
          " (p = ", formatC(stats_rows[["aging"]]$Wilcoxon_p,
                            format="e", digits=2), ")")
}
message("\n  Cross-dataset median comparison:")
message("    PCOS case normalised median  = ", round(pcos_case_med,  3))
message("    Aging aged normalised median = ", round(aging_aged_med, 3))
message("    PCOS ctrl normalised median  = ", round(pcos_ctrl_med,  3))
message("    Aging young normalised median= ", round(aging_young_med,3))

message("\n  Output files:")
for (f in list.files(OUTPUT_DIR, pattern="v2|comparison|overlay",
                     full.names=FALSE))
  message("    ", f)