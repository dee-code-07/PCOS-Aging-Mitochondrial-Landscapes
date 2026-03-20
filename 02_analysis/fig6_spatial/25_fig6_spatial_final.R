# spatial transcriptomics module scoring — complete script
# pcos (gse296728) and aging (gse188257, ya_1/ya_2=young, ya_3/ya_4=aged)

# packages
library(Seurat)
library(UCell)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(writexl)
library(scales)

# output directories
outdir <- "E:/Documents/mini_project/analysis/19_spatial_FINAL"
dir.create(file.path(outdir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables"),  recursive = TRUE, showWarnings = FALSE)

# gene sets
to_mouse <- function(genes) {
  paste0(toupper(substr(genes, 1, 1)),
         tolower(substr(genes, 2, nchar(genes))))
}

gene_sets <- list(
  Mitochondrial = to_mouse(read.csv(
    "E:/Documents/mini_project/scrna/resources/modules/mitochondrial_module.csv")$gene),
  Oxidative     = to_mouse(read.csv(
    "E:/Documents/mini_project/scrna/resources/modules/oxidative_stress_module.csv")$gene),
  Senescence    = to_mouse(read.csv(
    "E:/Documents/mini_project/scrna/resources/modules/senescence_SASP_module.csv")$gene),
  Inflammation  = to_mouse(read.csv(
    "E:/Documents/mini_project/scrna/resources/modules/inflammation_module.csv")$gene)
)

# sample paths
lt_base <- "E:/Documents/mini_project/spatial/output/06_label_transfer_FIXED"

# pcos: 3 controls + 3 pcos
pcos_samples <- list(
  list(rds = file.path(lt_base, "pcos/c1/c1_with_celltypes.rds"),
       label = "Control 1", condition = "Control"),
  list(rds = file.path(lt_base, "pcos/c2/c2_with_celltypes.rds"),
       label = "Control 2", condition = "Control"),
  list(rds = file.path(lt_base, "pcos/c3/c3_with_celltypes.rds"),
       label = "Control 3", condition = "Control"),
  list(rds = file.path(lt_base, "pcos/p1/p1_with_celltypes.rds"),
       label = "PCOS 1", condition = "PCOS"),
  list(rds = file.path(lt_base, "pcos/p2/p2_with_celltypes.rds"),
       label = "PCOS 2", condition = "PCOS"),
  list(rds = file.path(lt_base, "pcos/p3/p3_with_celltypes.rds"),
       label = "PCOS 3", condition = "PCOS")
)

# aging: ya_1 + ya_2 = young, ya_3 + ya_4 = aged (russ et al. 2022)
aging_samples <- list(
  list(rds = file.path(lt_base, "aging/ya_1/ya_1_with_celltypes.rds"),
       label = "Young 1", condition = "Young"),
  list(rds = file.path(lt_base, "aging/ya_2/ya_2_with_celltypes.rds"),
       label = "Young 2", condition = "Young"),
  list(rds = file.path(lt_base, "aging/ya_3/ya_3_with_celltypes.rds"),
       label = "Aged 3",  condition = "Aged"),
  list(rds = file.path(lt_base, "aging/ya_4/ya_4_with_celltypes.rds"),
       label = "Aged 4",  condition = "Aged")
)

# function: score one sample
score_sample <- function(samp_info, gene_sets) {
  
  if (!file.exists(samp_info$rds)) {
    message("MISSING: ", samp_info$rds); return(NULL)
  }
  
  so <- tryCatch(readRDS(samp_info$rds), error = function(e) {
    message("Load error: ", e$message); return(NULL)
  })
  if (is.null(so)) return(NULL)
  
  DefaultAssay(so) <- "RNA"
  obj_genes <- rownames(so)
  
  gs_filtered <- lapply(gene_sets, function(g) g[g %in% obj_genes])
  for (nm in names(gs_filtered)) {
    cat(samp_info$label, nm, ":", length(gs_filtered[[nm]]), "genes\n")
  }
  
  so <- tryCatch(
    AddModuleScore_UCell(so, features = gs_filtered, name = "_UCell"),
    error = function(e) { message("UCell error: ", e$message); return(NULL) }
  )
  if (is.null(so)) return(NULL)
  
  meta      <- so@meta.data
  ucell_cols <- grep("_UCell", colnames(meta), value = TRUE)
  
# pixel coordinates
  if ("pxl_col" %in% colnames(meta) && "pxl_row" %in% colnames(meta)) {
    meta$x <- as.numeric(meta$pxl_col)
    meta$y <- as.numeric(meta$pxl_row)
  } else if ("imagecol" %in% colnames(meta) && "imagerow" %in% colnames(meta)) {
    meta$x <- as.numeric(meta$imagecol)
    meta$y <- as.numeric(meta$imagerow)
  } else {
    stop("Cannot find pixel coordinates for: ", samp_info$label)
  }
  
  meta$sample_label <- samp_info$label
  meta$condition    <- samp_info$condition
  meta$barcode      <- rownames(meta)
  
  return(meta[, c("barcode", "x", "y", "sample_label", "condition", ucell_cols)])
}

# function: spatial map
plot_spatial_module <- function(all_meta, module_col, module_name,
                                dataset_label, ncols = 3) {
  df       <- all_meta %>% filter(!is.na(.data[[module_col]]))
  df$y_plot <- -df$y
  q_hi     <- quantile(df[[module_col]], 0.99, na.rm = TRUE)
  q_lo     <- min(df[[module_col]], na.rm = TRUE)
  
  plots <- list()
  for (sl in unique(df$sample_label)) {
    sub <- df[df$sample_label == sl, ]
    p <- ggplot(sub, aes(x = x, y = y_plot, colour = .data[[module_col]])) +
      geom_point(size = 0.6, alpha = 0.9) +
      scale_colour_gradientn(
        colours = c("#0D0887","#7E03A8","#CC4678","#F89441","#F0F921"),
        limits  = c(q_lo, q_hi), name = "UCell\nScore", oob = squish
      ) +
      coord_equal() + theme_void() + ggtitle(sl) +
      theme(plot.title      = element_text(face = "bold", size = 9, hjust = 0.5),
            legend.position = "right",
            legend.key.height = unit(0.8, "cm"),
            legend.title    = element_text(size = 8),
            legend.text     = element_text(size = 7))
    plots[[sl]] <- p
  }
  
  nrow_val <- ceiling(length(plots) / ncols)
  wrap_plots(plots, ncol = ncols, nrow = nrow_val) +
    plot_annotation(
      title    = paste0(module_name, " Score — ", dataset_label, " Spatial Dataset"),
      subtitle = paste0("Mus musculus ovary (Visium) | UCell scoring | ",
                        length(plots), " samples"),
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(colour = "grey50", size = 9)
      )
    )
}

# function: violin plot
theme_pub <- function() {
  theme_classic(base_size = 12) +
    theme(legend.position   = "bottom",
          legend.title      = element_blank(),
          strip.background  = element_rect(fill = "white", colour = "grey70"),
          strip.text        = element_text(face = "bold", size = 11),
          axis.text         = element_text(colour = "black"),
          plot.title        = element_text(face = "bold", size = 13),
          plot.subtitle     = element_text(colour = "grey50", size = 9))
}

make_violin <- function(meta, cond_levels, cond_labels,
                        fill_vals, title_str, ucell_cols) {
  df <- meta %>%
    select(condition, all_of(ucell_cols)) %>%
    pivot_longer(-condition, names_to = "Module", values_to = "Score") %>%
    mutate(Module    = gsub("_UCell.*", "", Module),
           condition = factor(condition,
                              levels = cond_levels, labels = cond_labels))
  
  ggplot(df, aes(x = condition, y = Score, fill = condition)) +
    geom_violin(alpha = 0.8, colour = "grey30",
                linewidth = 0.3, scale = "width") +
    geom_boxplot(width = 0.12, outlier.size = 0.2,
                 fill = "white", colour = "black", linewidth = 0.35) +
    scale_fill_manual(values = fill_vals) +
    facet_wrap(~ Module, scales = "free_y", nrow = 2) +
    labs(title    = title_str,
         subtitle = "Distribution across all tissue spots (Visium)",
         y = "UCell Score", x = NULL) +
    theme_pub()
}

# function: run wilcoxon stats
run_stats <- function(meta, g1_cond, g2_cond, comp_label, ucell_cols) {
  lapply(ucell_cols, function(col) {
    g1 <- meta[[col]][meta$condition == g1_cond]
    g2 <- meta[[col]][meta$condition == g2_cond]
    g1 <- g1[!is.na(g1)]; g2 <- g2[!is.na(g2)]
    if (length(g1) < 3 || length(g2) < 3) return(NULL)
    w  <- wilcox.test(g1, g2)
    d  <- (mean(g1) - mean(g2)) / sd(c(g1, g2))
    data.frame(
      comparison  = comp_label,
      module      = gsub("_UCell.*", "", col),
      n_group1    = length(g1),
      n_group2    = length(g2),
      mean_group1 = round(mean(g1), 5),
      mean_group2 = round(mean(g2), 5),
      cohens_d    = round(d, 4),
      p_value     = w$p.value
    )
  }) %>% bind_rows()
}

# process samples

message("\n=== Processing PCOS samples ===")
pcos_meta_list <- lapply(pcos_samples, score_sample, gene_sets = gene_sets)
pcos_meta_list <- pcos_meta_list[!sapply(pcos_meta_list, is.null)]
pcos_meta      <- bind_rows(pcos_meta_list)
cat("PCOS spots scored:", nrow(pcos_meta), "\n")

message("\n=== Processing Aging samples ===")
aging_meta_list <- lapply(aging_samples, score_sample, gene_sets = gene_sets)
aging_meta_list <- aging_meta_list[!sapply(aging_meta_list, is.null)]
aging_meta      <- bind_rows(aging_meta_list)
cat("Aging spots scored:", nrow(aging_meta), "\n")

# confirm ucell column names
ucell_cols <- grep("_UCell", colnames(pcos_meta), value = TRUE)
cat("\nUCell columns detected:", paste(ucell_cols, collapse = ", "), "\n")

cat("\nPCOS mean scores:\n")
print(colMeans(pcos_meta[, ucell_cols], na.rm = TRUE))
cat("\nAging mean scores:\n")
print(colMeans(aging_meta[, ucell_cols], na.rm = TRUE))

# statistics

message("\n=== Running statistics ===")

# pcos: pcos vs control
pcos_stats        <- run_stats(pcos_meta, "PCOS", "Control",
                               "PCOS_Case_vs_Control", ucell_cols)
pcos_stats$FDR    <- p.adjust(pcos_stats$p_value, "BH")
pcos_stats$sig    <- ifelse(pcos_stats$FDR < 0.001, "***",
                            ifelse(pcos_stats$FDR < 0.01,  "**",
                                   ifelse(pcos_stats$FDR < 0.05,  "*", "ns")))

cat("\nPCOS stats (PCOS vs Control):\n")
print(pcos_stats[, c("module","mean_group1","mean_group2","cohens_d","FDR","sig")])

# aging: aged vs young
aging_stats        <- run_stats(aging_meta, "Aged", "Young",
                                "Aged_vs_Young", ucell_cols)
aging_stats$FDR    <- p.adjust(aging_stats$p_value, "BH")
aging_stats$sig    <- ifelse(aging_stats$FDR < 0.001, "***",
                             ifelse(aging_stats$FDR < 0.01,  "**",
                                    ifelse(aging_stats$FDR < 0.05,  "*", "ns")))

cat("\nAging stats (Aged vs Young):\n")
print(aging_stats[, c("module","mean_group1","mean_group2","cohens_d","FDR","sig")])

# save table s11 — both sheets
write_xlsx(
  list("PCOS_module_stats"  = pcos_stats,
       "Aging_module_stats" = aging_stats),
  file.path(outdir, "tables", "Table_S11_spatial_module_stats_FINAL.xlsx")
)
message("Table S11 saved (both sheets).")

# spatial maps

message("\n=== Generating spatial maps ===")

mod_names <- gsub("_UCell.*", "", ucell_cols)

for (i in seq_along(ucell_cols)) {
  col  <- ucell_cols[i]
  name <- mod_names[i]
  
# pcos map (3x2 grid)
  tryCatch({
    p <- plot_spatial_module(pcos_meta, col, name, "PCOS", ncols = 3)
    fname <- file.path(outdir, "figures",
                       paste0("spatial_", name, "_PCOS_FINAL.tiff"))
    tiff(fname, width = 15, height = 9, units = "in",
         res = 600, compression = "lzw")
    print(p); dev.off()
    message(name, " PCOS map saved.")
  }, error = function(e) message("ERROR PCOS ", name, ": ", e$message))
  
# aging map (2x2 grid)
  tryCatch({
    p2 <- plot_spatial_module(aging_meta, col, name, "Aging", ncols = 2)
    fname2 <- file.path(outdir, "figures",
                        paste0("spatial_", name, "_Aging_FINAL.tiff"))
    tiff(fname2, width = 12, height = 11, units = "in",
         res = 600, compression = "lzw")
    print(p2); dev.off()
    message(name, " Aging map saved.")
  }, error = function(e) message("ERROR Aging ", name, ": ", e$message))
}

# violin plots

message("\n=== Generating violin plots ===")

# pcos violin: pcos vs control
p_pcos <- make_violin(
  pcos_meta,
  c("PCOS", "Control"), c("PCOS", "Control"),
  c("PCOS" = "#E07B39", "Control" = "#E8C06A"),
  "Spatial Module Scores: PCOS vs Control",
  ucell_cols
)
tiff(file.path(outdir, "figures", "spatial_violin_PCOS_FINAL.tiff"),
     width = 9, height = 7, units = "in", res = 600, compression = "lzw")
print(p_pcos); dev.off()
message("PCOS violin saved.")

# aging violin: aged vs young
p_aging <- make_violin(
  aging_meta,
  c("Aged", "Young"), c("Aged", "Young"),
  c("Aged" = "#2166AC", "Young" = "#92C5DE"),
  "Spatial Module Scores: Aged vs Young",
  ucell_cols
)
tiff(file.path(outdir, "figures", "spatial_violin_Aging_FINAL.tiff"),
     width = 9, height = 7, units = "in", res = 600, compression = "lzw")
print(p_aging); dev.off()
message("Aging violin saved.")

# done

message("\n=== ALL DONE ===")
cat("Figures saved to:", file.path(outdir, "figures"), "\n")
cat("Tables saved to:",  file.path(outdir, "tables"),  "\n")
cat("\nFiles generated:\n")
print(list.files(file.path(outdir, "figures")))
print(list.files(file.path(outdir, "tables")))