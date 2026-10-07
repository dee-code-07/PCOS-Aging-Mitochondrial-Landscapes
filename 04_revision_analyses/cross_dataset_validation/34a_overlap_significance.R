#!/usr/bin/env Rscript
# Local run of 22a

base <- "/mnt/e/mini_project/scrna/output/08_de"
f_pcos  <- file.path(base, "pcos/DE_pcos_case_vs_ctrl_universe.csv")
f_aging <- file.path(base, "aging/DE_aging_aged_vs_young_universe.csv")
out_dir <- "/mnt/e/mini_project/analysis/22_permutation_null/v3"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(42); B <- 10000
PADJ <- 0.05; LFC <- 0.25                                      # same as manuscript

rd <- function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  g <- if ("gene" %in% names(d)) d$gene else if (names(d)[1] %in% c("", "X")) d[[1]] else rownames(d)
  d$gene <- g; d[!duplicated(d$gene), ]
}
p <- rd(f_pcos); a <- rd(f_aging)
cat(sprintf("Rows (genes tested): PCOS %d, Aging %d\n", nrow(p), nrow(a)))

status <- function(d) ifelse(d$p_val_adj < PADJ & d$avg_log2FC >  LFC,  1L,
                      ifelse(d$p_val_adj < PADJ & d$avg_log2FC < -LFC, -1L, 0L))
p$s <- status(p); a$s <- status(a)

# Universe = genes tested in BOTH analyses (NOT the whole genome)
U <- intersect(p$gene, a$gene); N <- length(U)
ps <- p$s[match(U, p$gene)]; as <- a$s[match(U, a$gene)]
cat(sprintf("Universe N = %d | PCOS DEGs %d (up %d, down %d) | Aging DEGs %d (up %d, down %d)\n",
            N, sum(ps!=0), sum(ps==1), sum(ps==-1), sum(as!=0), sum(as==1), sum(as==-1)))

obs_up   <- sum(ps ==  1 & as ==  1)
obs_down <- sum(ps == -1 & as == -1)
obs      <- obs_up + obs_down
cat(sprintf("Observed same-direction overlap = %d (up %d, down %d); manuscript = 83 (39/44)\n",
            obs, obs_up, obs_down))
if (obs != 83) warning("Observed overlap != 83: thresholds/tables differ from the ones used for the manuscript. Resolve before using these numbers.")

hyp <- function(o, K, k) phyper(o - 1, K, N - K, k, lower.tail = FALSE)
res <- data.frame(
  test = c("UP", "DOWN"),
  observed = c(obs_up, obs_down),
  expected = c(sum(ps==1)*sum(as==1)/N, sum(ps==-1)*sum(as==-1)/N),
  p_hypergeom = c(hyp(obs_up, sum(ps==1), sum(as==1)), hyp(obs_down, sum(ps==-1), sum(as==-1))))

# Combined same-direction overlap: gene-level permutation (keeps DEG counts and directions)
null <- replicate(B, { s <- sample(as); sum(ps != 0 & ps == s) })
stopifnot(sd(null) > 0, !all(null == obs))            # sanity checks
emp_p <- (sum(null >= obs) + 1) / (B + 1)
res <- rbind(res, data.frame(test = "UP+DOWN (permutation)", observed = obs,
                             expected = mean(null), p_hypergeom = emp_p))
res$z <- c(NA, NA, (obs - mean(null)) / sd(null))
print(res, digits = 4)
cat(sprintf("Null: mean %.1f, SD %.2f, range %d-%d\n", mean(null), sd(null), min(null), max(null)))
