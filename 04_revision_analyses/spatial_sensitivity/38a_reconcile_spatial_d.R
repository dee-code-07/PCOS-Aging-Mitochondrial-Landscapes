#!/usr/bin/env Rscript
# =============================================================================
# 26b: Reconcile Script 19 vs Script 26 spatial d values
# 
# Goal: Load the EXACT same objects and gene sets as Script 19, compute
#       per-section means, and calculate d under BOTH formulas (sd(c()) and 
#       pooled SD) to resolve the discrepancy.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(UCell)
  library(dplyr)
})

set.seed(42)

outdir <- "/mnt/e/mini_project/analysis/26_spatial_sensitivity"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ── Gene sets: EXACT same as Script 19 ──────────────────────────────────────
to_mouse <- function(genes) {
  paste0(toupper(substr(genes, 1, 1)),
         tolower(substr(genes, 2, nchar(genes))))
}

mod_base <- "/mnt/e/mini_project/scrna/resources/modules"
gene_sets <- list(
  Mitochondrial = to_mouse(read.csv(
    file.path(mod_base, "mitochondrial_module.csv"))$gene),
  Oxidative     = to_mouse(read.csv(
    file.path(mod_base, "oxidative_stress_module.csv"))$gene),
  Senescence    = to_mouse(read.csv(
    file.path(mod_base, "senescence_SASP_module.csv"))$gene),
  Inflammation  = to_mouse(read.csv(
    file.path(mod_base, "inflammation_module.csv"))$gene)
)

# Print gene set sizes for audit
for (nm in names(gene_sets)) {
  cat(nm, ":", length(gene_sets[[nm]]), "genes\n")
  # Show first few mt genes if Mitochondrial
  if (nm == "Mitochondrial") {
    mt_genes <- grep("^[Mm]t", gene_sets[[nm]], value = TRUE)
    cat("  mt-prefixed genes:", paste(head(mt_genes, 10), collapse=", "), "\n")
  }
}

# ── Load samples: EXACT same as Script 19 ────────────────────────────────────
lt_base <- "/mnt/e/mini_project/spatial/output/06_label_transfer_FIXED"

aging_samples <- list(
  list(rds = file.path(lt_base, "aging/ya_1/ya_1_with_celltypes.rds"),
       label = "Young 1", section = "ya_1", condition = "Young"),
  list(rds = file.path(lt_base, "aging/ya_2/ya_2_with_celltypes.rds"),
       label = "Young 2", section = "ya_2", condition = "Young"),
  list(rds = file.path(lt_base, "aging/ya_3/ya_3_with_celltypes.rds"),
       label = "Aged 3", section = "ya_3", condition = "Aged"),
  list(rds = file.path(lt_base, "aging/ya_4/ya_4_with_celltypes.rds"),
       label = "Aged 4", section = "ya_4", condition = "Aged")
)

# Score each sample SEPARATELY (same approach as Script 19)
all_meta <- list()
for (samp in aging_samples) {
  cat("\nLoading:", samp$label, "\n")
  so <- readRDS(samp$rds)
  DefaultAssay(so) <- "RNA"
  
  obj_genes <- rownames(so)
  gs_filtered <- lapply(gene_sets, function(g) g[g %in% obj_genes])
  for (nm in names(gs_filtered)) {
    cat("  ", nm, ":", length(gs_filtered[[nm]]), "/", length(gene_sets[[nm]]), "genes matched\n")
  }
  
  so <- AddModuleScore_UCell(so, features = gs_filtered, name = "_UCell")
  
  meta <- so@meta.data
  ucell_cols <- grep("_UCell", colnames(meta), value = TRUE)
  
  meta$section <- samp$section
  meta$condition <- samp$condition
  meta$label <- samp$label
  
  all_meta[[samp$section]] <- meta[, c("section", "condition", "label", ucell_cols)]
  cat("  n_spots:", nrow(meta), "\n")
}

aging_meta <- bind_rows(all_meta)
cat("\nTotal aging spots:", nrow(aging_meta), "\n")

# ── Per-section summary ──────────────────────────────────────────────────────
ucell_cols <- grep("_UCell", colnames(aging_meta), value = TRUE)

cat("\n========== PER-SECTION MEANS ==========\n")
section_summary <- aging_meta %>%
  group_by(section, condition) %>%
  summarise(
    n_spots = n(),
    across(all_of(ucell_cols), mean, .names = "mean_{.col}"),
    .groups = "drop"
  )
print(as.data.frame(section_summary))

# Save per-section means
write.csv(section_summary, file.path(outdir, "per_section_means.csv"), row.names = FALSE)

# ── Cohen's d under BOTH formulas ────────────────────────────────────────────
cat("\n========== COHEN'S D: BOTH FORMULAS ==========\n")

d_total_sd <- function(g1, g2) {
  (mean(g1) - mean(g2)) / sd(c(g1, g2))
}

d_pooled_sd <- function(g1, g2) {
  n1 <- length(g1); n2 <- length(g2)
  pooled <- sqrt(((n1-1)*var(g1) + (n2-1)*var(g2)) / (n1+n2-2))
  (mean(g1) - mean(g2)) / pooled
}

# Combined aged vs young (same as Script 19)
results <- list()

for (col in ucell_cols) {
  g1 <- aging_meta[[col]][aging_meta$condition == "Aged"]
  g2 <- aging_meta[[col]][aging_meta$condition == "Young"]
  
  results[[length(results)+1]] <- data.frame(
    comparison = "Combined_Aged_vs_Young",
    module = gsub("_UCell.*", "", col),
    n_aged = length(g1),
    n_young = length(g2),
    mean_aged = round(mean(g1), 5),
    mean_young = round(mean(g2), 5),
    d_total_sd = round(d_total_sd(g1, g2), 4),
    d_pooled_sd = round(d_pooled_sd(g1, g2), 4),
    stringsAsFactors = FALSE
  )
}

# Per-section: ya_3 only vs Young, ya_4 only vs Young
for (aged_sec in c("ya_3", "ya_4")) {
  for (col in ucell_cols) {
    g1 <- aging_meta[[col]][aging_meta$section == aged_sec]
    g2 <- aging_meta[[col]][aging_meta$condition == "Young"]
    
    results[[length(results)+1]] <- data.frame(
      comparison = paste0(aged_sec, "_vs_Young"),
      module = gsub("_UCell.*", "", col),
      n_aged = length(g1),
      n_young = length(g2),
      mean_aged = round(mean(g1), 5),
      mean_young = round(mean(g2), 5),
      d_total_sd = round(d_total_sd(g1, g2), 4),
      d_pooled_sd = round(d_pooled_sd(g1, g2), 4),
      stringsAsFactors = FALSE
    )
  }
}

results_df <- bind_rows(results)

cat("\nFull results:\n")
print(results_df)

# Comparison with published values
cat("\n========== VERIFICATION ==========\n")
cat("Manuscript (Script 19):  Senescence d = +0.444, Mito d = -0.311\n")
sen_row <- results_df[results_df$comparison == "Combined_Aged_vs_Young" & 
                       results_df$module == "Senescence", ]
mito_row <- results_df[results_df$comparison == "Combined_Aged_vs_Young" & 
                        results_df$module == "Mitochondrial", ]
cat("This script (d_total_sd):", 
    "  Sen =", sen_row$d_total_sd, 
    "  Mito =", mito_row$d_total_sd, "\n")
cat("This script (d_pooled)::", 
    "  Sen =", sen_row$d_pooled_sd, 
    "  Mito =", mito_row$d_pooled_sd, "\n")
cat("Script 26 (d_pooled):  Sen = +0.464, Mito = -0.397\n")

# Save complete results
write.csv(results_df, file.path(outdir, "reconciled_d_values.csv"), row.names = FALSE)

# ── Section-pairwise: does each aged section beat each young section? ─────────
cat("\n========== SECTION-PAIRWISE COMPARISON ==========\n")
cat("Do BOTH aged sections score above BOTH young sections?\n\n")

for (col in ucell_cols) {
  mod_name <- gsub("_UCell.*", "", col)
  cat(mod_name, ":\n")
  for (sec in c("ya_1", "ya_2", "ya_3", "ya_4")) {
    vals <- aging_meta[[col]][aging_meta$section == sec]
    cat(sprintf("  %s (n=%d): mean=%.5f, median=%.5f, sd=%.5f\n",
                sec, length(vals), mean(vals), median(vals), sd(vals)))
  }
  # Check ordering
  means <- sapply(c("ya_1", "ya_2", "ya_3", "ya_4"), function(s) {
    mean(aging_meta[[col]][aging_meta$section == s])
  })
  both_aged_above <- means["ya_3"] > means["ya_1"] && 
                     means["ya_3"] > means["ya_2"] &&
                     means["ya_4"] > means["ya_1"] && 
                     means["ya_4"] > means["ya_2"]
  cat(sprintf("  All 4 pairwise: aged > young? %s\n\n",
              ifelse(both_aged_above, "YES ✓", "NO ✗")))
}

cat("\n========== DONE ==========\n")
