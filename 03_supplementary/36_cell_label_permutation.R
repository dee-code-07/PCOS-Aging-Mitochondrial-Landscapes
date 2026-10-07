#!/usr/bin/env Rscript
# 22b_cell_label_permutation.R  -- supplementary check only (run 22a first).
# Differences from the old script: explicit cells.1/cells.2 (no Idents reliance),
# same-direction statistic (as the 83 genes), per-permutation diagnostics,
# hard sanity aborts, checkpoints. SMOKE TEST FIRST: Rscript 22b... 3
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })
args <- commandArgs(TRUE); N_PERM <- if (length(args)) as.integer(args[1]) else 100
set.seed(42); options(future.globals.maxSize = Inf); future::plan("sequential")
out_dir <- "/home/deekshah/mini_project/analysis/22_permutation_null/v3"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

load_obj <- function(f, lab1, pat1, lab2, pat2) {
  o <- readRDS(f); DefaultAssay(o) <- "RNA"
  if (inherits(o[["RNA"]], "Assay5")) o[["RNA"]] <- JoinLayers(o[["RNA"]])
  o$condition <- ifelse(grepl(pat1, o$sample_id), lab1, ifelse(grepl(pat2, o$sample_id), lab2, NA))
  o <- subset(o, cells = colnames(o)[!is.na(o$condition)])
  keep <- unlist(lapply(c(lab1, lab2), function(l) {
    cc <- colnames(o)[o$condition == l]; if (length(cc) > 5000) sample(cc, 5000) else cc }))
  subset(o, cells = keep)
}
pcos  <- load_obj("/home/deekshah/mini_project/scrna/output/07_annotation/pcos/pcos_annotated.rds",
                  "PCOS", "case", "Control", "ctrl|control")
aging <- load_obj("/home/deekshah/mini_project/scrna/output/07_annotation/aging/aging_annotated.rds",
                  "Aged", "_a_|aged", "Young", "_y_|young")
cat("Cells per condition:\n"); print(table(pcos$condition)); print(table(aging$condition))

de <- function(o, c1, c2) {
  d <- FindMarkers(o, cells.1 = c1, cells.2 = c2, assay = "RNA", test.use = "wilcox",
                   min.pct = 0.1, verbose = FALSE)
  d$gene <- rownames(d)
  d$s <- ifelse(d$p_val_adj < 0.05 & d$avg_log2FC > 0.25, 1L,
          ifelse(d$p_val_adj < 0.05 & d$avg_log2FC < -0.25, -1L, 0L)); d
}
overlap <- function(p, a) { m <- inner_join(p[, c("gene","s")], a[, c("gene","s")], by = "gene")
                            sum(m$s.x != 0 & m$s.x == m$s.y) }
cells_of <- function(o, lab, v) colnames(o)[v == lab]

# Observed
p0 <- de(pcos, cells_of(pcos, "PCOS", pcos$condition), cells_of(pcos, "Control", pcos$condition))
a0 <- de(aging, cells_of(aging, "Aged", aging$condition), cells_of(aging, "Young", aging$condition))
obs <- overlap(p0, a0)
cat(sprintf("Observed same-direction overlap (subsampled): %d | DEGs PCOS %d, Aging %d\n",
            obs, sum(p0$s != 0), sum(a0$s != 0)))

null <- numeric(N_PERM); np <- numeric(N_PERM); na <- numeric(N_PERM)
for (i in seq_len(N_PERM)) {
  sp <- sample(pcos$condition); sa <- sample(aging$condition)   # shuffled label vectors
  agree <- mean(sp == pcos$condition)                           # should be ~0.5
  if (agree > 0.9) stop("Shuffle did not change labels (agreement ", agree, ")")
  p <- de(pcos,  cells_of(pcos,  "PCOS", sp), cells_of(pcos,  "Control", sp))
  a <- de(aging, cells_of(aging, "Aged", sa), cells_of(aging, "Young",   sa))
  null[i] <- overlap(p, a); np[i] <- sum(p$s != 0); na[i] <- sum(a$s != 0)
  cat(sprintf("perm %d: label agreement %.2f | DEGs PCOS %d Aging %d | overlap %d\n",
              i, agree, np[i], na[i], null[i]))
  if (i == 5 && all(null[1:5] == obs)) stop("First 5 null overlaps equal the observed value: label use is broken.")
  closeAllConnections(); gc(verbose = FALSE)
  if (i %% 10 == 0) write.csv(data.frame(perm = 1:i, overlap = null[1:i], n_pcos = np[1:i], n_aging = na[1:i]),
                              file.path(out_dir, "cell_perm_checkpoint.csv"), row.names = FALSE)
}
if (sd(null) == 0) stop("Null SD = 0: permutation invalid.")
if (any(null == obs)) warning("At least one permutation equals the observed value.")
emp_p <- (sum(null >= obs) + 1) / (N_PERM + 1)
write.csv(data.frame(perm = 1:N_PERM, overlap = null, n_pcos = np, n_aging = na),
          file.path(out_dir, "cell_permutation_null.csv"), row.names = FALSE)
write.csv(data.frame(observed = obs, mean_null = mean(null), sd_null = sd(null),
                     empirical_p = emp_p, z = (obs - mean(null)) / sd(null)),
          file.path(out_dir, "cell_permutation_results.csv"), row.names = FALSE)
cat(sprintf("Obs %d | null mean %.2f SD %.2f | empirical p %.4f\n", obs, mean(null), sd(null), emp_p))
