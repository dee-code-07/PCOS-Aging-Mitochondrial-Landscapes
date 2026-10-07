suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
  library(ggplot2)
  library(patchwork)
})

in_csv <- "/mnt/e/mini_project/analysis/ai_revision_tasks/pseudotime_data.csv"
out_dir <- "/mnt/e/mini_project/revision_final/s19_v2"
old_xlsx <- "/mnt/e/mini_project/revisions/Supplementary_Tables/Supplementary_Table_S19.xlsx"

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

# Load data
pt <- read.csv(in_csv)
pt$condition_group <- ifelse(grepl("ageseq_a", pt$sample_id), "Aging_aged",
                           ifelse(grepl("ageseq_y", pt$sample_id), "Aging_young",
                           ifelse(grepl("case", pt$sample_id), "PCOS_case",
                           ifelse(grepl("ctrl", pt$sample_id), "PCOS_control", "Other"))))

# Compute Medians and Pairwise
get_summary <- function(col) {
  pt %>% filter(!is.na(!!sym(col))) %>%
    group_by(condition_group) %>%
    summarize(median_pseudotime = median(!!sym(col)), n_cells = n()) %>%
    as.data.frame()
}

get_pairwise <- function(col) {
  sub_pt <- pt[!is.na(pt[[col]]), ]
  pw <- pairwise.wilcox.test(sub_pt[[col]], sub_pt$condition_group, p.adjust.method="BH")
  pw_df <- as.data.frame(as.table(pw$p.value))
  colnames(pw_df) <- c("Group1", "Group2", "p.adj")
  pw_df <- pw_df %>% filter(!is.na(p.adj))
  pw_df
}

harm_meds <- get_summary("pseudotime_harmony")
harm_pw <- get_pairwise("pseudotime_harmony")

rpca_meds <- get_summary("pseudotime_rpca")
rpca_pw <- get_pairwise("pseudotime_rpca")

# Spearman rho
valid_cells <- pt[!is.na(pt$pseudotime_rpca) & !is.na(pt$pseudotime_harmony), ]
rho <- cor(valid_cells$pseudotime_rpca, valid_cells$pseudotime_harmony, method="spearman")
n_rho <- nrow(valid_cells)
write.csv(data.frame(spearman_rho=rho, n_cells=n_rho), file.path(out_dir, "spearman_correlation.csv"), row.names=FALSE)

# Load old sheets
old_rpca_meds <- read.xlsx(old_xlsx, sheet="C_RPCA_medians")
old_rpca_pw <- read.xlsx(old_xlsx, sheet="D_RPCA_pairwise")
old_harm_meds <- read.xlsx(old_xlsx, sheet="A_Harmony_medians")
old_harm_pw <- read.xlsx(old_xlsx, sheet="B_Harmony_pairwise")

# Check RPCA
check_meds <- all.equal(rpca_meds, old_rpca_meds)
check_pw <- all.equal(rpca_pw, old_rpca_pw)
cat("RPCA medians match existing: ", isTRUE(check_meds), "\n")
cat("RPCA pairwise match existing: ", isTRUE(check_pw), "\n")

# Write new XLSX
wb <- loadWorkbook(old_xlsx)
removeWorksheet(wb, "A_Harmony_medians")
removeWorksheet(wb, "B_Harmony_pairwise")
addWorksheet(wb, "A_Harmony_medians")
addWorksheet(wb, "B_Harmony_pairwise")
writeData(wb, "A_Harmony_medians", harm_meds)
writeData(wb, "B_Harmony_pairwise", harm_pw)
# Reorder sheets to keep the README first and group them reasonably
worksheetOrder(wb) <- c(1, 4, 5, 2, 3) # README, A_H, B_H, C_R, D_R
saveWorkbook(wb, file.path(out_dir, "Supplementary_Table_S19_v2.xlsx"), overwrite=TRUE)

# Plot Figure S19
plot_pt <- function(df, col, title, meds) {
  ggplot(df[!is.na(df[[col]]), ], aes(x=condition_group, y=!!sym(col), fill=condition_group)) +
    geom_violin(alpha=0.6, trim=TRUE) +
    geom_point(data=meds, aes(x=condition_group, y=median_pseudotime, fill=condition_group), shape=23, size=3, color="black", fill="black") +
    theme_bw() + labs(title=title, y="Pseudotime", x="") +
    theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1))
}
p1 <- plot_pt(pt, "pseudotime_harmony", "Harmony Pseudotime", harm_meds)
p2 <- plot_pt(pt, "pseudotime_rpca", "RPCA Pseudotime", rpca_meds)
fig <- p1 + p2
ggsave(file.path(out_dir, "Figure_S19_v2.png"), fig, width=10, height=5, dpi=300)
ggsave(file.path(out_dir, "Figure_S19_v2.pdf"), fig, width=10, height=5)

cat("\n=== Harmony Medians (New vs Old) ===\n")
comp_meds <- merge(harm_meds, old_harm_meds, by="condition_group", suffixes=c("_New", "_Old"))
print(comp_meds)

cat("\n=== Harmony Pairwise (New vs Old) ===\n")
comp_pw <- merge(harm_pw, old_harm_pw, by=c("Group1", "Group2"), suffixes=c("_New", "_Old"), all=TRUE)
print(comp_pw)

cat("\n=== Spearman Rho ===\n")
cat("Rho:", rho, "N:", n_rho, "\n")
