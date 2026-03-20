#!/usr/bin/env Rscript
# 04_module_statistics_fixed.r
# what was broken in the enhanced version:
# - used orig.ident for condition assignment — wrong, all cells got same label
# - rstatix::wilcox_test + cohens_d silently failed -> all stats na
# - lme4 mixed models wrote empty files
# - violin significance annotations failed silently
# - heatmap was from old pcos-vs-aging script, not within-dataset comparisons
# fixes:
# - uses sample_id (established pattern from pseudotime/deg scripts)
# - pcos:  grepl("ctrl|control") → control  |  grepl("case") → pcos
# - aging: grepl("_y_|young")    → young    |  grepl("_a_|aged|old") → aged
# - base r wilcox.test() — no rstatix dependency
# - manual cohen's d pooled-sd formula
# - cell type column: final_celltype with seurat_clusters fallback
# - significance annotations via annotate() — no rstatix needed
# input:  analysis/03_modules/pcos_with_modules.rds
# analysis/03_modules/aging_with_modules.rds
# output: analysis/04_module_statistics/

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(patchwork)
})

# paths
PROJECT_ROOT <- "E:/Documents/mini_project"
PCOS_INPUT   <- file.path(PROJECT_ROOT, "analysis/03_modules/pcos_with_modules.rds")
AGING_INPUT  <- file.path(PROJECT_ROOT, "analysis/03_modules/aging_with_modules.rds")
OUTPUT_DIR   <- file.path(PROJECT_ROOT, "analysis/04_module_statistics")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(PCOS_INPUT))  stop("Missing: ", PCOS_INPUT)
if (!file.exists(AGING_INPUT)) stop("Missing: ", AGING_INPUT)

cat("Loading data...\n")
pcos  <- readRDS(PCOS_INPUT)
aging <- readRDS(AGING_INPUT)

# module detection
keywords <- c("Mitochondrial", "Senescence", "Oxidative", "Inflammation")
modules <- map_chr(keywords, function(k) {
  m <- grep(paste0(k, ".*1$"), colnames(pcos@meta.data), value = TRUE,
            ignore.case = TRUE)
  if (length(m) == 0) stop("Module score not found for keyword: ", k)
  m[1]
})
names(modules) <- keywords
cat("Modules detected:", paste(modules, collapse = ", "), "\n")

# condition assignment — uses sample_id (established project convention)
# pattern confirmed across: pseudotime script (l3443-3458), deg script (l3739),
# and celltype deconvolution scripts throughout analysis.txt

# diagnostic: confirm sample_id exists
if (!"sample_id" %in% colnames(pcos@meta.data))
  stop("'sample_id' not found in PCOS metadata. Available cols: ",
       paste(colnames(pcos@meta.data), collapse=", "))
if (!"sample_id" %in% colnames(aging@meta.data))
  stop("'sample_id' not found in Aging metadata.")

cat("\nUnique PCOS sample_ids:  ", paste(sort(unique(pcos$sample_id)),  collapse=", "), "\n")
cat("Unique Aging sample_ids: ", paste(sort(unique(aging$sample_id)), collapse=", "), "\n")

# pcos: ctrl/control → control  |  case → pcos
pcos_meta <- pcos@meta.data %>%
  mutate(
    condition = case_when(
      grepl("ctrl|control", sample_id, ignore.case = TRUE) ~ "Control",
      grepl("case",         sample_id, ignore.case = TRUE) ~ "PCOS",
      TRUE ~ "Unknown"
    ),
    celltype = if ("final_celltype" %in% colnames(.)) final_celltype
    else as.character(seurat_clusters)
  ) %>%
  filter(condition %in% c("PCOS", "Control"))

# aging: _y_/young → young  |  _a_/aged/old → aged
aging_meta <- aging@meta.data %>%
  mutate(
    condition = case_when(
      grepl("_y_|young",      sample_id, ignore.case = TRUE) ~ "Young",
      grepl("_a_|aged|old",   sample_id, ignore.case = TRUE) ~ "Aged",
      TRUE ~ "Unknown"
    ),
    celltype = if ("final_celltype" %in% colnames(.)) final_celltype
    else as.character(seurat_clusters)
  ) %>%
  filter(condition %in% c("Aged", "Young"))

# diagnostic
cat("\nPCOS condition table:\n");  print(table(pcos_meta$condition))
cat("\nAging condition table:\n"); print(table(aging_meta$condition))

if (length(unique(pcos_meta$condition))  < 2)
  stop("PCOS: only one condition found.\n",
       "sample_id values seen: ", paste(unique(pcos_meta$sample_id), collapse=", "))
if (length(unique(aging_meta$condition)) < 2)
  stop("Aging: only one condition found.\n",
       "sample_id values seen: ", paste(unique(aging_meta$sample_id), collapse=", "))

# statistical functions

# cohen's d — pooled sd formula
cohens_d_manual <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 2 || length(y) < 2) return(NA_real_)
  pooled_sd <- sqrt(((length(x)-1)*var(x) + (length(y)-1)*var(y)) /
                      (length(x) + length(y) - 2))
  if (pooled_sd == 0) return(NA_real_)
  (mean(x) - mean(y)) / pooled_sd
}

# wilcoxon + cohen's d for one module, one data slice
run_test <- function(df, mod, label_a, label_b) {
  x <- df[[mod]][df$condition == label_a]
  y <- df[[mod]][df$condition == label_b]
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 3 || length(y) < 3) {
    return(data.frame(n_a=length(x), n_b=length(y),
                      mean_a=NA, mean_b=NA, W=NA,
                      p_value=NA, p_adj=NA,
                      cohens_d=NA, direction=NA,
                      stringsAsFactors=FALSE))
  }
  wt  <- wilcox.test(x, y, exact = FALSE)
  cd  <- cohens_d_manual(x, y)
  data.frame(
    n_a       = length(x),
    n_b       = length(y),
    mean_a    = round(mean(x), 5),
    mean_b    = round(mean(y), 5),
    W         = as.numeric(wt$statistic),
    p_value   = wt$p.value,
    p_adj     = NA_real_,
    cohens_d  = round(cd, 4),
    direction = ifelse(cd > 0,
                       paste0(label_a, "_higher"),
                       paste0(label_b, "_higher")),
    stringsAsFactors = FALSE
  )
}

# overall stats
calc_overall_stats <- function(meta, label_a, label_b, dataset_label) {
  res <- map_dfr(names(modules), function(k) {
    mod <- modules[[k]]
    r   <- run_test(meta, mod, label_a, label_b)
    r$module  <- mod
    r$stratum <- "Overall"
    r$dataset <- dataset_label
    r
  })
  res$p_adj <- p.adjust(res$p_value, method = "BH")
  res %>% select(dataset, module, stratum, n_a, n_b,
                 mean_a, mean_b, W, p_value, p_adj, cohens_d, direction)
}

# cell-type stratified stats
calc_celltype_stats <- function(meta, label_a, label_b, dataset_label) {
  cts <- unique(meta$celltype)
  res <- map_dfr(cts, function(ct) {
    sub <- meta %>% filter(celltype == ct)
    map_dfr(names(modules), function(k) {
      mod <- modules[[k]]
      r   <- run_test(sub, mod, label_a, label_b)
      r$module  <- mod
      r$stratum <- ct
      r$dataset <- dataset_label
      r
    })
  })
  res <- res %>% filter(!is.na(p_value))
  if (nrow(res) > 0) res$p_adj <- p.adjust(res$p_value, method = "BH")
  res %>% select(dataset, module, stratum, n_a, n_b,
                 mean_a, mean_b, W, p_value, p_adj, cohens_d, direction)
}

# run stats
cat("\nRunning PCOS statistics (PCOS vs Control)...\n")
pcos_overall  <- calc_overall_stats(pcos_meta,  "PCOS",  "Control", "PCOS")
pcos_celltype <- calc_celltype_stats(pcos_meta, "PCOS",  "Control", "PCOS")

cat("Running Aging statistics (Aged vs Young)...\n")
aging_overall  <- calc_overall_stats(aging_meta,  "Aged",  "Young", "Aging")
aging_celltype <- calc_celltype_stats(aging_meta, "Aged",  "Young", "Aging")

# write stats
write_csv(pcos_overall,   file.path(OUTPUT_DIR, "PCOS_stats_overall.csv"))
write_csv(pcos_celltype,  file.path(OUTPUT_DIR, "PCOS_stats_celltype.csv"))
write_csv(aging_overall,  file.path(OUTPUT_DIR, "Aging_stats_overall.csv"))
write_csv(aging_celltype, file.path(OUTPUT_DIR, "Aging_stats_celltype.csv"))

# pseudobulk per sample
pseudobulk <- bind_rows(
  pcos_meta  %>% group_by(sample_id, condition) %>%
    summarise(across(all_of(unname(modules)), mean, na.rm=TRUE), .groups="drop") %>%
    mutate(dataset = "PCOS"),
  aging_meta %>% group_by(sample_id, condition) %>%
    summarise(across(all_of(unname(modules)), mean, na.rm=TRUE), .groups="drop") %>%
    mutate(dataset = "Aging")
)
write_csv(pseudobulk, file.path(OUTPUT_DIR, "pseudobulk_per_sample.csv"))

cat("\nPCOS overall results:\n")
print(pcos_overall[, c("module","p_value","p_adj","cohens_d","direction")])
cat("\nAging overall results:\n")
print(aging_overall[, c("module","p_value","p_adj","cohens_d","direction")])

# theme
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text    = element_text(colour = "black"),
      axis.title.x = element_blank(),
      plot.title   = element_text(face = "bold", hjust = 0.5,
                                  size = base_size + 1),
      legend.position = "none"
    )
}

pval_label <- function(p) {
  if (is.na(p))  return("n.s.")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("n.s.")
}

# violin plots
plot_violin <- function(meta, mod, title, fill_vals, stats_row, prefix) {
  pd <- meta %>%
    select(condition, score = all_of(mod)) %>%
    filter(!is.na(score))
  
  y_bar   <- quantile(pd$score, 0.98, na.rm=TRUE) * 1.10
  y_label <- y_bar * 1.06
  conds   <- sort(unique(pd$condition))
  p_lab   <- pval_label(stats_row$p_adj)
  cd_lab  <- if (!is.na(stats_row$cohens_d))
    sprintf("d = %.2f", stats_row$cohens_d) else ""
  
  p <- ggplot(pd, aes(x = condition, y = score, fill = condition)) +
    geom_violin(trim = FALSE, alpha = 0.85, width = 0.8) +
    geom_boxplot(width = 0.12, fill = "white",
                 outlier.shape = NA, linewidth = 0.5) +
    scale_fill_manual(values = fill_vals) +
    labs(title = title, y = "Module score") +
    theme_pub()
  
  if (!is.na(stats_row$p_value) && length(conds) == 2) {
    p <- p +
      annotate("segment", x=1, xend=2, y=y_bar, yend=y_bar, linewidth=0.5) +
      annotate("text", x=1.5, y=y_label,
               label = paste0(p_lab, "\n", cd_lab),
               size=3.2, hjust=0.5)
  }
  
  fname <- file.path(OUTPUT_DIR,
                     paste0(prefix, "_", gsub("_Score1","",mod), "_Violin.tiff"))
  tiff(fname, width=4, height=5, units="in", res=600, compression="lzw")
  print(p); dev.off()
  message("  Saved: ", basename(fname))
}

cat("\nGenerating violin plots...\n")
pcos_fills  <- c(PCOS    = "#D55E00", Control = "#0072B2")
aging_fills <- c(Aged    = "#CC79A7", Young   = "#009E73")

for (mod in unname(modules)) {
  sr_p <- pcos_overall  %>% filter(module == mod)
  sr_a <- aging_overall %>% filter(module == mod)
  if (nrow(sr_p) == 0) sr_p <- data.frame(p_value=NA, p_adj=NA, cohens_d=NA)
  if (nrow(sr_a) == 0) sr_a <- data.frame(p_value=NA, p_adj=NA, cohens_d=NA)
  
  plot_violin(pcos_meta,  mod,
              title     = paste("PCOS —",  gsub("_Score1","",mod)),
              fill_vals = pcos_fills,
              stats_row = sr_p[1,], prefix = "PCOS")
  
  plot_violin(aging_meta, mod,
              title     = paste("Aging —", gsub("_Score1","",mod)),
              fill_vals = aging_fills,
              stats_row = sr_a[1,], prefix = "Aging")
}

# cell-type × module heatmap
cat("\nGenerating cell-type heatmap...\n")

heat_data <- bind_rows(pcos_celltype, aging_celltype) %>%
  filter(!is.na(cohens_d)) %>%
  mutate(
    module_label = gsub("_Score1", "", module),
    sig_label    = ifelse(!is.na(p_adj) & p_adj < 0.05, "*", "")
  )

if (nrow(heat_data) > 0) {
  ct_order <- heat_data %>%
    group_by(stratum) %>%
    summarise(mean_abs = mean(abs(cohens_d), na.rm=TRUE)) %>%
    arrange(desc(mean_abs)) %>% pull(stratum)
  
  heat_data$stratum      <- factor(heat_data$stratum, levels = rev(ct_order))
  heat_data$module_label <- factor(heat_data$module_label)
  heat_data$dataset      <- factor(heat_data$dataset, levels = c("PCOS","Aging"))
  
  lim <- ceiling(max(abs(heat_data$cohens_d), na.rm=TRUE) * 10) / 10
  
  ph <- ggplot(heat_data,
               aes(x = module_label, y = stratum, fill = cohens_d)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_text(aes(label = sig_label), size = 4.5,
              colour = "black", vjust = 0.75) +
    scale_fill_gradient2(
      low = "#4682B4", mid = "white", high = "#B22222",
      midpoint = 0, limits = c(-lim, lim),
      name = "Cohen's d\n(positive = disease/aged higher)"
    ) +
    facet_wrap(~ dataset, scales = "free_x") +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x      = element_text(angle=35, hjust=1, size=10, face="bold"),
      axis.text.y      = element_text(size=9),
      axis.title       = element_blank(),
      strip.text       = element_text(face="bold", size=12),
      strip.background = element_rect(fill="grey92", colour=NA),
      legend.title     = element_text(size=9),
      plot.title       = element_text(face="bold", hjust=0.5, size=13),
      plot.subtitle    = element_text(hjust=0.5, colour="grey40", size=10)
    ) +
    labs(
      title    = "Module Score Effect Sizes by Cell Type",
      subtitle = "Cohen's d (disease vs control)  |  * FDR < 0.05"
    )
  
  n_ct <- length(unique(heat_data$stratum))
  tiff(file.path(OUTPUT_DIR, "Module_celltype_heatmap_600dpi.tiff"),
       width=10, height=max(5, 0.35*n_ct + 2.5),
       units="in", res=600, compression="lzw")
  print(ph); dev.off()
  message("  Saved: Module_celltype_heatmap_600dpi.tiff")
} else {
  message("  WARNING: No cell-type stats for heatmap — check celltype column exists.")
}

# summary
cat("\n", strrep("=",60), "\n  OUTPUT SUMMARY\n", strrep("=",60), "\n")
cat("  PCOS_stats_overall.csv   —", nrow(pcos_overall),   "rows\n")
cat("  PCOS_stats_celltype.csv  —", nrow(pcos_celltype),  "rows\n")
cat("  Aging_stats_overall.csv  —", nrow(aging_overall),  "rows\n")
cat("  Aging_stats_celltype.csv —", nrow(aging_celltype), "rows\n")
cat("  pseudobulk_per_sample.csv\n")
cat("  8 violin TIFFs\n")
cat("  Module_celltype_heatmap_600dpi.tiff\n")
cat(strrep("=",60), "\nDone.\n")