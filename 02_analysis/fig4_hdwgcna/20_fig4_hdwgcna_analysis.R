# hdwgcna analysis: single-cell co-expression network analysis
# project: single-cell and spatial landscapes of mitochondrial dysfunction
# in pcos and ovarian aging
# author:  [your name]
# date:    2026
# biological rationale:
# pcos and ovarian aging share overlapping metabolic, oxidative, and
# inflammatory perturbations. mitochondria are a central hub connecting
# these pathways. hdwgcna (morabito et al., cell reports methods 2023)
# enables the identification of co-expressed gene modules in individual
# cell types from sparse single-cell data via metacell aggregation.
# we run hdwgcna separately in granulosa, stromal, and immune compartments
# to reveal cell-type-specific mitochondrial dysfunction and senescence
# programs shared between pcos and aged ovaries.
# references:
# 1. morabito s et al. hdwgcna identifies co-expression networks in
# high-dimensional transcriptomics data. cell rep methods. 2023.
# https://doi.org/10.1016/j.crmeth.2023.100498
# 2. langfelder p & horvath s. wgcna: an r package for weighted
# correlation network analysis. bmc bioinformatics. 2008.
# https://doi.org/10.1186/1471-2105-9-559
# 3. hao y et al. integrated analysis of multimodal single-cell data.
# cell. 2021. https://doi.org/10.1016/j.cell.2021.04.048

# package installation (run once if needed)
# uncomment and run if packages are not yet installed:
# install.packages(c("seurat", "tidyverse", "cowplot", "patchwork",
# "ggplot2", "viridis", "rcolorbrewer", "scales"))
# install.packages("biocmanager")
# biocmanager::install(c("wgcna", "ucell", "geneoverlap"))
# install.packages("devtools")
# devtools::install_github("smorabit/hdwgcna", ref = "dev")
# install.packages("harmony")
# biocmanager::install("enrichr")  # for enrichment analysis

# load libraries
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
  library(harmony)
  library(UCell)
  library(future)
  library(parallel)
})

# global ggplot theme for publication quality
theme_set(theme_cowplot(font_size = 12))

# random seed for full reproducibility
set.seed(42)

# completely disable all parallel backends
plan(sequential)

options(future.fork.enable = FALSE)
options(mc.cores = 1)
options(parallelly.fork.enable = FALSE)

# reset wgcna threading to safe value
allowWGCNAThreads(nThreads = 2)

# define paths
base_dir    <- "E:/Documents/mini_project"
scrna_dir   <- file.path(base_dir, "scrna/output")
out_dir     <- file.path(base_dir, "hdwgcna")
plot_dir    <- file.path(out_dir, "plots")
rds_dir     <- file.path(out_dir, "rds")

# create output directories
for (d in c(out_dir, plot_dir, rds_dir,
            file.path(plot_dir, "granulosa"),
            file.path(plot_dir, "stromal"),
            file.path(plot_dir, "immune"))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# input rds paths
aging_rds <- file.path(scrna_dir, "07_annotation/aging/aging_annotated.rds")
pcos_rds  <- file.path(scrna_dir, "07_annotation/pcos/pcos_annotated.rds")

cat("✓ Paths defined.\n")
cat("  Output directory:", out_dir, "\n\n")

# load and merge seurat objects
# we merge the pcos and aging datasets into a single seurat object.
# this is essential for integrated co-expression analysis that can reveal
# shared and divergent gene modules across conditions.
# we keep orig.ident as "aging"/"pcos" and use sample_id for
# biological replicate tracking (critical for metacell construction).

cat("Loading aging dataset...\n")
aging_obj <- readRDS(aging_rds)
aging_obj <- SeuratObject::UpdateSeuratObject(aging_obj)

cat("Loading PCOS dataset...\n")
pcos_obj <- readRDS(pcos_rds)
pcos_obj <- SeuratObject::UpdateSeuratObject(pcos_obj)

# standardize: ensure rna assay is default before merge
DefaultAssay(aging_obj) <- "RNA"
DefaultAssay(pcos_obj)  <- "RNA"

# critical: remove sct and integrated assays before merging.
# seurat v5 cannot rbind sct residual matrices of different dimensions
# (different numbers of variable genes per dataset). we only need the
# raw rna counts + normalized data for hdwgcna downstream.
aging_obj[["SCT"]]        <- NULL
aging_obj[["integrated"]] <- NULL
pcos_obj[["SCT"]]         <- NULL
pcos_obj[["integrated"]]  <- NULL
gc()
cat("✓ SCT/integrated assays removed. Proceeding with RNA only.\n")

# cell type harmonization
# both datasets use 'final_celltype'. we create a 'broad_celltype' column
# that collapses subtypes into three compartments for hdwgcna grouping.
# subtype information is preserved in 'final_celltype' for downstream use.

# granulosa subtypes (all follicular stages)
granulosa_labels <- c("granulosa", "granulosa_antral",
                      "granulosa_luteal", "granulosa_preantral")

# stromal subtypes (mesenchymal + vascular)
stromal_labels <- c("stroma", "stroma_fibroblast", "stroma_pericyte",
                    "theca", "endothelium", "lymphatic_endo",
                    "epithelium", "thelium")

# immune subtypes
immune_labels <- c("b_cells", "t_cells", "nk", "phagocytes",
                   "dendritic", "immune")

# function to map final_celltype → broad_celltype
assign_broad <- function(seurat_obj) {
  ct <- seurat_obj$final_celltype
  broad <- case_when(
    ct %in% granulosa_labels ~ "Granulosa",
    ct %in% stromal_labels   ~ "Stromal",
    ct %in% immune_labels    ~ "Immune",
    ct == "oocyte"           ~ "Oocyte",
    TRUE                     ~ "Other"
  )
  seurat_obj$broad_celltype <- broad
  return(seurat_obj)
}

aging_obj <- assign_broad(aging_obj)
pcos_obj  <- assign_broad(pcos_obj)

# add condition labels derived from sample_id
# aging: ageseq_a_* = aged, ageseq_y_* = youngcontrol
# pcos:  pcos_case_* = pcos_case, pcos_ctrl_* = pcos_control
assign_condition <- function(seurat_obj) {
  sid <- seurat_obj$sample_id
  cond <- case_when(
    grepl("ageseq_a", sid) ~ "Aged",
    grepl("ageseq_y", sid) ~ "YoungControl",
    grepl("pcos_case", sid) ~ "PCOS_Case",
    grepl("pcos_ctrl", sid) ~ "PCOS_Control",
    TRUE ~ "Unknown"
  )
  seurat_obj$condition <- cond
  return(seurat_obj)
}

aging_obj <- assign_condition(aging_obj)
pcos_obj  <- assign_condition(pcos_obj)

# verify broad_celltype assignment
cat("\n=== AGING: broad_celltype distribution ===\n")
print(table(aging_obj$broad_celltype))
cat("\n=== PCOS: broad_celltype distribution ===\n")
print(table(pcos_obj$broad_celltype))

# merge
# merge both objects. we do not re-integrate here; hdwgcna handles
# batch correction via metacell grouping by sample_id.
cat("\nMerging datasets...\n")
seurat_merged <- merge(
  x = aging_obj,
  y = pcos_obj,
  add.cell.ids = c("aging", "pcos"),
  project = "PCOS_Aging_hdWGCNA",
  merge.data = TRUE
)

# clean up individual objects to free ram
rm(aging_obj, pcos_obj)
gc()

# ensure rna assay active and join layers (seurat v5 requirement)
DefaultAssay(seurat_merged) <- "RNA"
seurat_merged <- JoinLayers(seurat_merged)

cat("✓ Merge complete.\n")
cat("  Total cells:", ncol(seurat_merged), "\n")
cat("  Conditions:", paste(unique(seurat_merged$condition), collapse = ", "), "\n")
cat("  Broad cell types:", paste(unique(seurat_merged$broad_celltype), collapse = ", "), "\n\n")

# re-process merged object
# hdwgcna requires normalizedata, findvariablefeatures, scaledata, pca, umap
# on the merged object so gene selection and metacell knn are valid.

cat("Re-processing merged object (Normalize → HVG → Scale → PCA → Harmony → UMAP)...\n")
cat("This may take 5-10 minutes on 16GB RAM...\n")

seurat_merged <- NormalizeData(seurat_merged, normalization.method = "LogNormalize",
                               scale.factor = 10000, verbose = FALSE)

seurat_merged <- FindVariableFeatures(seurat_merged, selection.method = "vst",
                                      nfeatures = 3000, verbose = FALSE)

# scale only variable features to save memory
seurat_merged <- ScaleData(seurat_merged,
                           features = VariableFeatures(seurat_merged),
                           verbose = FALSE)

seurat_merged <- RunPCA(seurat_merged,
                        features = VariableFeatures(seurat_merged),
                        npcs = 40, verbose = FALSE)

# harmony batch correction across sample_id (biological replicates)
seurat_merged <- RunHarmony(seurat_merged,
                            group.by.vars = "sample_id",
                            reduction = "pca",
                            reduction.save = "harmony",
                            verbose = FALSE)

seurat_merged <- RunUMAP(seurat_merged,
                         reduction = "harmony",
                         dims = 1:30,
                         reduction.name = "umap",
                         verbose = FALSE)

cat("✓ Re-processing complete.\n\n")

# baseline qc / overview plots
# plot 1: umap by broad cell type
# plot 2: umap by condition
# plot 3: umap by final_celltype (all subtypes)
# plot 4: cell count bar plot per condition x broad_celltype

cat("Generating overview plots...\n")

# color palettes
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

# plot 1: umap – broad cell type
p_umap_broad <- DimPlot(seurat_merged,
                        group.by = "broad_celltype",
                        cols = broad_colors,
                        pt.size = 0.3,
                        label = TRUE,
                        label.size = 4,
                        repel = TRUE) +
  umap_theme() +
  ggtitle("Merged Ovarian Atlas\nBroad Cell Type") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 4)))

# plot 2: umap – condition
p_umap_cond <- DimPlot(seurat_merged,
                       group.by = "condition",
                       cols = condition_colors,
                       pt.size = 0.3,
                       label = FALSE) +
  umap_theme() +
  ggtitle("Merged Ovarian Atlas\nCondition") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        legend.position = "right") +
  guides(color = guide_legend(override.aes = list(size = 4)))

# plot 3: umap – final cell type (all subtypes)
p_umap_final <- DimPlot(seurat_merged,
                        group.by = "final_celltype",
                        pt.size = 0.2,
                        label = TRUE,
                        label.size = 3,
                        repel = TRUE) +
  umap_theme() +
  ggtitle("Merged Ovarian Atlas\nAll Cell Subtypes") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        legend.position = "right") +
  NoLegend()

# combine and save
p_overview <- (p_umap_broad | p_umap_cond) / p_umap_final +
  plot_annotation(
    title = "PCOS and Ovarian Aging – Integrated Single-Cell Atlas",
    subtitle = paste0("n = ", format(ncol(seurat_merged), big.mark=","), " cells | PCOS (case/ctrl) + Aging (young/aged)"),
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                  plot.subtitle = element_text(size = 11, hjust = 0.5))
  )

ggsave(file.path(plot_dir, "01_UMAP_merged_overview.pdf"),
       p_overview, width = 16, height = 14, dpi = 300)
ggsave(file.path(plot_dir, "01_UMAP_merged_overview.png"),
       p_overview, width = 16, height = 14, dpi = 300)

# plot 4: cell composition bar chart
cell_comp <- seurat_merged@meta.data %>%
  filter(broad_celltype %in% c("Granulosa", "Stromal", "Immune")) %>%
  group_by(condition, broad_celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(pct = n / sum(n) * 100)

p_comp <- ggplot(cell_comp, aes(x = condition, y = pct, fill = broad_celltype)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  scale_fill_manual(values = broad_colors, name = "Cell Type") +
  scale_x_discrete(limits = c("YoungControl", "Aged", "PCOS_Control", "PCOS_Case")) +
  labs(title = "Cell Type Composition Across Conditions",
       subtitle = "Percentage of Granulosa, Stromal and Immune cells",
       x = "Condition", y = "Cell Proportion (%)") +
  theme_cowplot(font_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

ggsave(file.path(plot_dir, "02_cell_composition_barplot.pdf"),
       p_comp, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, "02_cell_composition_barplot.png"),
       p_comp, width = 8, height = 6, dpi = 300)

cat("✓ Overview plots saved.\n\n")

# save merged object checkpoint
cat("Saving merged Seurat object...\n")
saveRDS(seurat_merged, file.path(rds_dir, "seurat_merged_preprocessed.rds"))
cat("✓ Saved to:", file.path(rds_dir, "seurat_merged_preprocessed.rds"), "\n\n")

# hdwgcna setup function
# we define a master function that runs the complete hdwgcna pipeline for
# a given broad cell type compartment. this ensures:
# - consistent methodology across granulosa, stromal, immune
# - all tutorial plots generated for each compartment
# - memory-efficient sequential execution

run_hdwgcna <- function(seurat_obj,
                        cell_type_name,    # "Granulosa", "Stromal", "Immune"
                        wgcna_name,        # short name for the experiment
                        plot_subdir,       # subfolder in plot_dir
                        k_value = 25,      # metacell k (tuned per compartment)
                        fraction = 0.05,   # min gene expression fraction
                        soft_power = NULL, # NULL = auto-select
                        min_cells = 50) {  # min cells per group for metacell
  
  cat("\n", strrep("═", 70), "\n")
  cat("  hdWGCNA:", cell_type_name, "Compartment\n")
  cat(strrep("═", 70), "\n\n")
  
  pdir <- file.path(plot_dir, plot_subdir)
  if (!dir.exists(pdir)) dir.create(pdir, recursive = TRUE)
  
# setupforwgcna
# select genes expressed in ≥ fraction of cells.
# 'fraction' method is preferred over 'variable' for scrna-seq because it
# captures lowly-expressed but biologically important genes (e.g., mt-genes)
  cat("[1/12] Setting up hdWGCNA experiment:", wgcna_name, "...\n")
  
  seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select  = "fraction",
    fraction     = fraction,
    wgcna_name   = wgcna_name
  )
  cat("  Genes selected:", length(GetWGCNAGenes(seurat_obj, wgcna_name)), "\n")
  
# construct metacells
# metacells aggregate k similar cells from the same biological sample.
# grouping by both broad_celltype and sample_id ensures:
# (a) metacells are cell-type-pure
# (b) metacells don't mix samples (critical for biological validity)
  cat("[2/12] Constructing metacells (k =", k_value, ")...\n")
  
  seurat_obj <- MetacellsByGroups(
    seurat_obj   = seurat_obj,
    group.by     = c("broad_celltype", "sample_id"),
    ident.group  = "broad_celltype",
    reduction    = "harmony",
    k            = k_value,
    max_shared   = round(k_value * 0.4),  # 40% max overlap
    min_cells    = min_cells,
    wgcna_name   = wgcna_name,
    verbose      = FALSE
  )
  
  seurat_obj <- NormalizeMetacells(seurat_obj, wgcna_name = wgcna_name)
  cat("  ✓ Metacells constructed.\n")
  
# optional: visualize metacell umap
  tryCatch({
    seurat_obj <- ScaleMetacells(seurat_obj,
                                 features = VariableFeatures(seurat_obj),
                                 wgcna_name = wgcna_name)
    seurat_obj <- RunPCAMetacells(seurat_obj,
                                  features = VariableFeatures(seurat_obj),
                                  wgcna_name = wgcna_name)
    seurat_obj <- RunHarmonyMetacells(seurat_obj,
                                      group.by.vars = "sample_id",
                                      wgcna_name = wgcna_name,
                                      verbose = FALSE)
    seurat_obj <- RunUMAPMetacells(seurat_obj,
                                   reduction = "harmony",
                                   dims = 1:15,
                                   wgcna_name = wgcna_name,
                                   verbose = FALSE)
    
    p_mc1 <- DimPlotMetacells(seurat_obj,
                              group.by = "broad_celltype",
                              wgcna_name = wgcna_name) +
      umap_theme() +
      scale_color_manual(values = broad_colors) +
      ggtitle(paste0(cell_type_name, " Metacells\nCell Type")) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    
    p_mc2 <- DimPlotMetacells(seurat_obj,
                              group.by = "condition",
                              wgcna_name = wgcna_name) +
      umap_theme() +
      scale_color_manual(values = condition_colors) +
      ggtitle(paste0(cell_type_name, " Metacells\nCondition")) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5))
    
    p_mc <- p_mc1 | p_mc2
    ggsave(file.path(pdir, paste0("01_metacell_UMAP_", wgcna_name, ".pdf")),
           p_mc, width = 12, height = 5, dpi = 300)
    ggsave(file.path(pdir, paste0("01_metacell_UMAP_", wgcna_name, ".png")),
           p_mc, width = 12, height = 5, dpi = 300)
    cat("  ✓ Metacell UMAP saved.\n")
  }, error = function(e) {
    cat("  ⚠ Metacell UMAP skipped:", conditionMessage(e), "\n")
  })
  
# setdatexpr
  cat("[3/12] Setting expression matrix for", cell_type_name, "...\n")
  seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name    = cell_type_name,
    group.by      = "broad_celltype",
    assay         = "RNA",
    layer         = "data",
    use_metacells = TRUE,
    wgcna_name    = wgcna_name
  )
  cat("  ✓ Expression matrix set.\n")
  
# soft power threshold
# we test powers 1:30. scale-free topology r² ≥ 0.80 is the cutoff.
# a signed network is preferred for biological data as it distinguishes
# positive and negative co-expression, relevant for oxidative stress
# signaling where some genes up-regulate as others down-regulate.
  cat("[4/12] Testing soft power thresholds...\n")
  seurat_obj <- TestSoftPowers(
    seurat_obj,
    networkType = "signed",
    wgcna_name  = wgcna_name
  )
  
# get power table and auto-select if not provided
  power_table <- GetPowerTable(seurat_obj, wgcna_name = wgcna_name)
  
  if (is.null(soft_power)) {
# select lowest power with sft.r.sq >= 0.80
    good_powers <- power_table %>%
      filter(SFT.R.sq >= 0.80, slope < 0) %>%
      arrange(Power)
    if (nrow(good_powers) > 0) {
      soft_power <- good_powers$Power[1]
    } else {
# fallback: select power at max r²
      soft_power <- power_table$Power[which.max(power_table$SFT.R.sq)]
      cat("  ⚠ R² < 0.80 at all powers. Using power at max R²:", soft_power, "\n")
    }
  }
  cat("  ✓ Soft power selected:", soft_power, "\n")
  
# plot soft power
  plot_list <- PlotSoftPowers(seurat_obj, wgcna_name = wgcna_name)
  p_sp <- wrap_plots(plot_list, ncol = 2) +
    plot_annotation(
      title = paste0(cell_type_name, " – Soft Power Threshold Selection"),
      subtitle = paste0("Selected power = ", soft_power,
                        " | Signed network | Scale-free topology R² ≥ 0.80"),
      theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
                    plot.subtitle = element_text(size = 10, hjust = 0.5))
    )
  ggsave(file.path(pdir, paste0("02_softpower_", wgcna_name, ".pdf")),
         p_sp, width = 12, height = 8, dpi = 300)
  ggsave(file.path(pdir, paste0("02_softpower_", wgcna_name, ".png")),
         p_sp, width = 12, height = 8, dpi = 300)
  cat("  ✓ Soft power plot saved.\n")
  
# construct co-expression network
# blockwiseconsensusmodules is used under the hood.
# minmodulesize = 50 ensures biologically meaningful modules.
# mergecutheight = 0.25 merges very similar modules.
  
# free memory before most ram-intensive step
  gc()
  cat("  ✓ Memory cleaned before network construction.\n")
  
  cat("[5/12] Constructing co-expression network...\n")
  seurat_obj <- ConstructNetwork(
    seurat_obj,
    soft_power      = soft_power,
    networkType     = "signed",
    minModuleSize   = 50,
    mergeCutHeight  = 0.25,
    deepSplit       = 4,
    detectCutHeight = 0.995,
    wgcna_name      = wgcna_name,
    overwrite_tom   = TRUE
  )
  cat("  ✓ Network constructed.\n")
  
# dendrogram plot
  tryCatch({
    p_dend <- PlotDendrogram(seurat_obj,
                             main = paste0(cell_type_name,
                                           " – Gene Co-expression Dendrogram\n",
                                           "(Colors = Module Assignment)"),
                             wgcna_name = wgcna_name)
    ggsave(file.path(pdir, paste0("03_dendrogram_", wgcna_name, ".pdf")),
           p_dend, width = 12, height = 5, dpi = 300)
    ggsave(file.path(pdir, paste0("03_dendrogram_", wgcna_name, ".png")),
           p_dend, width = 12, height = 5, dpi = 300)
    cat("  ✓ Dendrogram saved.\n")
  }, error = function(e) {
    cat("  ⚠ Dendrogram plot error:", conditionMessage(e), "\n")
  })
  
# module eigengenes
# module eigengenes (mes) summarize each module's expression as a single
# value per cell — the first principal component of all module genes.
# harmonize_me = true corrects for sample-level batch effects.
  cat("[6/12] Computing module eigengenes...\n")
  seurat_obj <- ModuleEigengenes(
    seurat_obj,
    group.by.vars  = "sample_id",
    wgcna_name     = wgcna_name,
    verbose        = FALSE
  )
  cat("  ✓ Module eigengenes computed.\n")
  
# module connectivity
# kme (module membership) = pearson correlation between a gene's expression
# and its module's eigengene. hub genes = highest kme.
# these are biologically the most representative genes for each module.
  cat("[7/12] Computing module connectivity (hub genes)...\n")
  seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by   = "broad_celltype",
    group_name = cell_type_name,
    wgcna_name = wgcna_name
  )
  cat("  ✓ Module connectivity computed.\n")
  
# rename modules with biologically interpretable prefixes
  seurat_obj <- ResetModuleNames(
    seurat_obj,
    new_name   = paste0(substr(cell_type_name, 1, 3), "-M"),
    wgcna_name = wgcna_name
  )
  
# visualize modules on umap
  cat("[8/12] Generating module UMAP visualizations...\n")
  
# get module eigengene names
  MEs <- GetMEs(seurat_obj, harmonized = TRUE, wgcna_name = wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name = wgcna_name)
  module_names <- levels(modules$module)
  module_names <- module_names[module_names != "grey"]
  
# add mes to metadata for plotting
  seurat_obj@meta.data <- cbind(seurat_obj@meta.data,
                                MEs[rownames(seurat_obj@meta.data), , drop = FALSE])
  
# plot top modules on umap (max 9 per figure for readability)
  if (length(module_names) > 0) {
    top_mods <- module_names[1:min(9, length(module_names))]
    tryCatch({
      p_feats <- ModuleFeaturePlot(
        seurat_obj,
        features       = top_mods,
        order          = TRUE,
        ucell          = TRUE,
        wgcna_name     = wgcna_name
      )
      p_feat_combined <- wrap_plots(p_feats, ncol = 3) +
        plot_annotation(
          title = paste0(cell_type_name, " – Module Eigengene UMAP"),
          subtitle = "Color intensity = module eigengene score per cell",
          theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
                        plot.subtitle = element_text(size = 10, hjust = 0.5))
        )
      ggsave(file.path(pdir, paste0("04_module_featureplot_", wgcna_name, ".pdf")),
             p_feat_combined, width = 15, height = 5 * ceiling(length(top_mods)/3), dpi = 300)
      ggsave(file.path(pdir, paste0("04_module_featureplot_", wgcna_name, ".png")),
             p_feat_combined, width = 15, height = 5 * ceiling(length(top_mods)/3), dpi = 300)
      cat("  ✓ Module feature UMAP saved.\n")
    }, error = function(e) {
      cat("  ⚠ ModuleFeaturePlot error:", conditionMessage(e), "\n")
    })
  }
  
# hub gene dot plot
  cat("[9/12] Generating hub gene dot plot...\n")
  tryCatch({
    p_hub <- HubGeneNetworkPlot(
      seurat_obj,
      n_hubs     = 3,
      n_other    = 5,
      edge_prop  = 0.75,
      mods       = module_names,
      wgcna_name = wgcna_name
    )
    ggsave(file.path(pdir, paste0("05_hub_gene_network_", wgcna_name, ".pdf")),
           p_hub, width = 14, height = 12, dpi = 300)
    ggsave(file.path(pdir, paste0("05_hub_gene_network_", wgcna_name, ".png")),
           p_hub, width = 14, height = 12, dpi = 300)
    cat("  ✓ Hub gene network plot saved.\n")
  }, error = function(e) {
    cat("  ⚠ HubGeneNetworkPlot error:", conditionMessage(e), "\n")
  })
  
# module eigengene dot plot (by condition)
  cat("[10/12] Plotting module eigengenes by condition...\n")
  tryCatch({
# subset to this cell type only
    cells_keep <- which(seurat_obj$broad_celltype == cell_type_name)
    seurat_sub <- seurat_obj[, cells_keep]
    
    p_me_dot <- DotPlot(
      seurat_sub,
      features = colnames(MEs),
      group.by = "condition",
      assay    = "RNA"
    ) +
      coord_flip() +
      scale_color_gradient2(low = "#2166AC", mid = "white", high = "#D6604D",
                            midpoint = 0, name = "Avg. Expression") +
      labs(title = paste0(cell_type_name, " – Module Eigengene Activity by Condition"),
           x = "Module", y = "Condition") +
      theme_cowplot(font_size = 11) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            plot.title = element_text(face = "bold", hjust = 0.5))
    
    ggsave(file.path(pdir, paste0("06_ME_dotplot_condition_", wgcna_name, ".pdf")),
           p_me_dot, width = 10, height = max(6, length(colnames(MEs)) * 0.4 + 2), dpi = 300)
    ggsave(file.path(pdir, paste0("06_ME_dotplot_condition_", wgcna_name, ".png")),
           p_me_dot, width = 10, height = max(6, length(colnames(MEs)) * 0.4 + 2), dpi = 300)
    cat("  ✓ ME dot plot by condition saved.\n")
    rm(seurat_sub)
  }, error = function(e) {
    cat("  ⚠ ME dot plot error:", conditionMessage(e), "\n")
  })
  
# module-trait correlation
# correlate module eigengenes with biological traits:
# condition (pcos vs aging), percent.mt (mito stress proxy)
  cat("[11/12] Module-trait correlation analysis...\n")
  tryCatch({
# encode traits as numeric
    seurat_obj$is_PCOS_Case    <- as.numeric(seurat_obj$condition == "PCOS_Case")
    seurat_obj$is_Aged         <- as.numeric(seurat_obj$condition == "Aged")
    seurat_obj$is_PCOS_Control <- as.numeric(seurat_obj$condition == "PCOS_Control")
    
    traits <- c("is_PCOS_Case", "is_Aged", "is_PCOS_Control", "percent.mt")
    
    seurat_obj <- ModuleTraitCorrelation(
      seurat_obj,
      traits     = traits,
      group.by   = "broad_celltype",
      wgcna_name = wgcna_name
    )
    
    p_mtc <- PlotModuleTraitCorrelation(
      seurat_obj,
      label        = "fdr",
      label_symbol = "stars",
      trait_names  = c("PCOS Case", "Aged", "PCOS Control", "Mito Stress (%mt)"),
      plot_max     = 0.4,
      combine      = TRUE,
      wgcna_name   = wgcna_name
    )
    ggsave(file.path(pdir, paste0("07_module_trait_correlation_", wgcna_name, ".pdf")),
           p_mtc, width = 10, height = max(8, length(module_names) * 0.4 + 3), dpi = 300)
    ggsave(file.path(pdir, paste0("07_module_trait_correlation_", wgcna_name, ".png")),
           p_mtc, width = 10, height = max(8, length(module_names) * 0.4 + 3), dpi = 300)
    cat("  ✓ Module-trait correlation plot saved.\n")
  }, error = function(e) {
    cat("  ⚠ Module-trait correlation error:", conditionMessage(e), "\n")
  })
  
# differential module eigengene (dme) analysis
# dme = statistical test for whether a module is significantly up- or
# down-regulated between two conditions. this is the core biological result:
# which modules are dysregulated in pcos_case vs pcos_control,
# and whether those same modules are also altered in aged vs youngcontrol.
  cat("[12/12] Differential Module Eigengene (DME) analysis...\n")
  tryCatch({
    cells_ct <- which(seurat_obj$broad_celltype == cell_type_name)
    
# pcos case vs control
    dme_pcos <- FindDMEs(
      seurat_obj,
      barcodes1  = names(cells_ct)[seurat_obj$condition[cells_ct] == "PCOS_Case"],
      barcodes2  = names(cells_ct)[seurat_obj$condition[cells_ct] == "PCOS_Control"],
      test.use   = "wilcox",
      wgcna_name = wgcna_name
    )
    dme_pcos$comparison <- "PCOS_Case_vs_Control"
    
# aged vs youngcontrol
    dme_aging <- FindDMEs(
      seurat_obj,
      barcodes1  = names(cells_ct)[seurat_obj$condition[cells_ct] == "Aged"],
      barcodes2  = names(cells_ct)[seurat_obj$condition[cells_ct] == "YoungControl"],
      test.use   = "wilcox",
      wgcna_name = wgcna_name
    )
    dme_aging$comparison <- "Aged_vs_YoungControl"
    
# save dme tables
    write.csv(dme_pcos,
              file.path(out_dir, paste0("DME_PCOS_", wgcna_name, ".csv")),
              row.names = FALSE)
    write.csv(dme_aging,
              file.path(out_dir, paste0("DME_Aging_", wgcna_name, ".csv")),
              row.names = FALSE)
    
# plot dme volcano-style
    dme_all <- rbind(dme_pcos, dme_aging)
    
    p_dme <- ggplot(dme_all, aes(x = avg_log2FC, y = -log10(p_val_adj),
                                 color = comparison, size = abs(avg_log2FC))) +
      geom_point(alpha = 0.85) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
      geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "grey40") +
      ggrepel::geom_text_repel(
        data = dme_all %>%
          filter(p_val_adj < 0.05) %>%
          group_by(comparison) %>%
          slice_max(abs(avg_log2FC), n = 5),
        aes(label = module),
        size = 3, max.overlaps = 20, show.legend = FALSE
      ) +
      scale_color_manual(values = c("PCOS_Case_vs_Control" = "#E63946",
                                    "Aged_vs_YoungControl"  = "#457B9D"),
                         name = "Comparison") +
      scale_size_continuous(range = c(2, 6), guide = "none") +
      labs(title = paste0(cell_type_name, " – Differential Module Eigengenes"),
           subtitle = "Significant modules (FDR < 0.05) shared between PCOS and Aging reveal convergent dysfunction",
           x = "Average log2 Fold Change",
           y = "-log10(Adjusted p-value)") +
      theme_cowplot(font_size = 12) +
      theme(plot.title    = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 9, hjust = 0.5),
            legend.position = "bottom")
    
    ggsave(file.path(pdir, paste0("08_DME_volcano_", wgcna_name, ".pdf")),
           p_dme, width = 10, height = 8, dpi = 300)
    ggsave(file.path(pdir, paste0("08_DME_volcano_", wgcna_name, ".png")),
           p_dme, width = 10, height = 8, dpi = 300)
    cat("  ✓ DME volcano plot saved.\n")
    
# dme heatmap: shared modules
    shared_sig <- intersect(
      dme_pcos %>% filter(p_val_adj < 0.05) %>% pull(module),
      dme_aging %>% filter(p_val_adj < 0.05) %>% pull(module)
    )
    if (length(shared_sig) > 0) {
      cat("  ★ SHARED SIGNIFICANT MODULES (PCOS ∩ Aging):", paste(shared_sig, collapse = ", "), "\n")
      write.csv(data.frame(shared_module = shared_sig,
                           cell_type = cell_type_name),
                file.path(out_dir, paste0("shared_modules_", wgcna_name, ".csv")),
                row.names = FALSE)
    }
    
    rm(dme_all, dme_pcos, dme_aging, cells_ct)
    gc()
  }, error = function(e) {
    cat("  ⚠ DME analysis error:", conditionMessage(e), "\n")
  })
  
# enrichment analysis
# go and kegg enrichment on each module's genes reveals biological functions.
# focus on mitochondrial, oxidative stress, and senescence terms.
  cat("  [Bonus] Running GO enrichment on modules...\n")
  tryCatch({
    seurat_obj <- RunEnrichr(
      seurat_obj,
      dbs        = c("GO_Biological_Process_2023",
                     "KEGG_2021_Human",
                     "MSigDB_Hallmark_2020",
                     "Reactome_2022"),
      max_genes  = 100,
      wgcna_name = wgcna_name
    )
    
    enrich_table <- GetEnrichrTable(seurat_obj, wgcna_name = wgcna_name)
    write.csv(enrich_table,
              file.path(out_dir, paste0("enrichment_", wgcna_name, ".csv")),
              row.names = FALSE)
    
# enrichment dot plot — top pathways per module
    p_enrich <- EnrichrDotPlot(
      seurat_obj,
      mods         = module_names,
      database     = "GO_Biological_Process_2023",
      n_terms      = 3,
      wgcna_name   = wgcna_name
    )
    ggsave(file.path(pdir, paste0("09_enrichment_dotplot_", wgcna_name, ".pdf")),
           p_enrich, width = 14, height = max(8, length(module_names) * 0.6 + 3), dpi = 300)
    ggsave(file.path(pdir, paste0("09_enrichment_dotplot_", wgcna_name, ".png")),
           p_enrich, width = 14, height = max(8, length(module_names) * 0.6 + 3), dpi = 300)
    cat("  ✓ Enrichment dot plot saved.\n")
  }, error = function(e) {
    cat("  ⚠ Enrichment error (network required):", conditionMessage(e), "\n")
    cat("    → Enrichment tables saved if available; plots skipped.\n")
  })
  
# module scores on umap (ucell)
  cat("  [Bonus] Scoring mitochondrial/senescence gene sets on cells...\n")
  tryCatch({
# core gene sets for biological story
    mito_genes <- c("MT-CO1","MT-CO2","MT-CO3","MT-ND1","MT-ND2","MT-ND4",
                    "MT-ND5","MT-ATP6","MT-CYB",
                    "TFAM","POLG","PINK1","PRKN","MFN1","MFN2","DNM1L",
                    "OPA1","FIS1","BNIP3","NIX")
    sasp_genes  <- c("IL6","IL8","IL1A","IL1B","CXCL1","CXCL2","CXCL8",
                     "MMP3","MMP9","MMP13","SERPINE1","IGFBP3","IGFBP7",
                     "CCL2","CCL20","TNF","VEGFA")
    senescence_genes <- c("CDKN1A","CDKN2A","TP53","RB1","LMNB1","HMGA1",
                          "HMGA2","DEC1","GLB1","PTPN11","MDM2")
    ros_genes   <- c("SOD1","SOD2","CAT","GPX1","GPX4","PRDX1","PRDX3",
                     "TXNRD2","NFE2L2","KEAP1","HMOX1","NQO1")
    
    gene_sets <- list(
      MitoStress  = mito_genes,
      SASP        = sasp_genes,
      Senescence  = senescence_genes,
      OxStress    = ros_genes
    )
    
    seurat_obj <- AddModuleScore_UCell(seurat_obj, features = gene_sets,
                                       name = NULL, ncores = 2)
    
    score_cols <- paste0(names(gene_sets), "_UCell")
    score_cols <- score_cols[score_cols %in% colnames(seurat_obj@meta.data)]
    
    if (length(score_cols) > 0) {
      p_scores <- lapply(score_cols, function(sc) {
        FeaturePlot(seurat_obj, features = sc,
                    pt.size = 0.2, order = TRUE) +
          scale_color_viridis(option = "magma", name = "Score") +
          umap_theme() +
          ggtitle(gsub("_UCell", "", sc)) +
          theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11))
      })
      
      p_scores_combined <- wrap_plots(p_scores, ncol = 2) +
        plot_annotation(
          title    = paste0(cell_type_name, " – Mitochondrial & Senescence Gene Set Scores"),
          subtitle = "UCell scoring | Mito stress, SASP, Senescence markers, Oxidative stress",
          theme    = theme(plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
                           plot.subtitle = element_text(size = 10, hjust = 0.5))
        )
      
      ggsave(file.path(pdir, paste0("10_geneset_scores_UMAP_", wgcna_name, ".pdf")),
             p_scores_combined, width = 12, height = 6 * ceiling(length(score_cols)/2), dpi = 300)
      ggsave(file.path(pdir, paste0("10_geneset_scores_UMAP_", wgcna_name, ".png")),
             p_scores_combined, width = 12, height = 6 * ceiling(length(score_cols)/2), dpi = 300)
      cat("  ✓ Gene set score UMAPs saved.\n")
      
# violin plot: scores per condition
      score_long <- seurat_obj@meta.data %>%
        filter(broad_celltype == cell_type_name) %>%
        dplyr::select(condition, all_of(score_cols)) %>%
        pivot_longer(cols = all_of(score_cols),
                     names_to = "GeneSet",
                     values_to = "Score") %>%
        mutate(GeneSet = gsub("_UCell", "", GeneSet),
               condition = factor(condition,
                                  levels = c("YoungControl","Aged","PCOS_Control","PCOS_Case")))
      
      p_vln <- ggplot(score_long, aes(x = condition, y = Score, fill = condition)) +
        geom_violin(trim = TRUE, scale = "width", alpha = 0.85) +
        geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.6) +
        facet_wrap(~GeneSet, scales = "free_y", ncol = 2) +
        scale_fill_manual(values = condition_colors) +
        stat_compare_means(comparisons = list(c("PCOS_Case","PCOS_Control"),
                                              c("Aged","YoungControl"),
                                              c("PCOS_Case","Aged")),
                           method = "wilcox.test",
                           label  = "p.signif",
                           hide.ns = TRUE) +
        labs(title    = paste0(cell_type_name, " – Gene Set Scores by Condition"),
             subtitle = "Violin + boxplot | Wilcoxon test (* p<0.05, ** p<0.01, *** p<0.001)",
             x = "Condition", y = "UCell Score") +
        theme_cowplot(font_size = 11) +
        theme(axis.text.x  = element_text(angle = 30, hjust = 1),
              strip.text   = element_text(face = "bold"),
              legend.position = "none",
              plot.title    = element_text(face = "bold", hjust = 0.5),
              plot.subtitle = element_text(size = 9, hjust = 0.5))
      
      ggsave(file.path(pdir, paste0("11_geneset_violin_condition_", wgcna_name, ".pdf")),
             p_vln, width = 10, height = 8, dpi = 300)
      ggsave(file.path(pdir, paste0("11_geneset_violin_condition_", wgcna_name, ".png")),
             p_vln, width = 10, height = 8, dpi = 300)
      cat("  ✓ Gene set violin plots saved.\n")
    }
  }, error = function(e) {
    cat("  ⚠ Gene set scoring error:", conditionMessage(e), "\n")
  })
  
# top hub gene heatmap
  cat("  [Bonus] Generating hub gene expression heatmap...\n")
  tryCatch({
    mods_df    <- GetModules(seurat_obj, wgcna_name = wgcna_name)
    hub_genes  <- mods_df %>%
      filter(module != "grey") %>%
      group_by(module) %>%
      slice_max(kME, n = 5) %>%
      ungroup() %>%
      arrange(module)
    
    cells_sub <- which(seurat_obj$broad_celltype == cell_type_name)
    seurat_sub <- seurat_obj[hub_genes$gene_name, cells_sub]
    
# order cells by condition
    cond_order <- c("YoungControl", "Aged", "PCOS_Control", "PCOS_Case")
    cell_order <- order(match(seurat_sub$condition, cond_order))
    seurat_sub <- seurat_sub[, cell_order]
    
    mat <- as.matrix(GetAssayData(seurat_sub, layer = "data"))
    mat <- t(scale(t(mat)))
    mat[mat > 2.5]  <- 2.5
    mat[mat < -2.5] <- -2.5
    
# annotation colors
    ann_col <- data.frame(
      Condition = seurat_sub$condition,
      row.names = colnames(seurat_sub)
    )
    ann_colors <- list(Condition = condition_colors[names(condition_colors) %in% unique(ann_col$Condition)])
    
# row annotation (module)
    mod_palette <- setNames(
      colorRampPalette(brewer.pal(9, "Set1"))(length(unique(hub_genes$module))),
      unique(hub_genes$module)
    )
    ann_row <- data.frame(Module = hub_genes$module, row.names = hub_genes$gene_name)
    ann_colors$Module <- mod_palette
    
    p_heat <- pheatmap::pheatmap(
      mat,
      annotation_col  = ann_col,
      annotation_row  = ann_row,
      annotation_colors = ann_colors,
      show_colnames   = FALSE,
      show_rownames   = TRUE,
      cluster_rows    = FALSE,
      cluster_cols    = FALSE,
      fontsize_row    = 7,
      color           = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      main            = paste0(cell_type_name, " – Top Hub Gene Expression\n(z-scored, ordered by condition)"),
      silent          = TRUE
    )
    ggsave(file.path(pdir, paste0("12_hub_gene_heatmap_", wgcna_name, ".pdf")),
           p_heat, width = 14, height = max(8, nrow(mat) * 0.25 + 3), dpi = 300)
    ggsave(file.path(pdir, paste0("12_hub_gene_heatmap_", wgcna_name, ".png")),
           p_heat, width = 14, height = max(8, nrow(mat) * 0.25 + 3), dpi = 300)
    cat("  ✓ Hub gene heatmap saved.\n")
    rm(seurat_sub, mat, hub_genes, mods_df, cells_sub)
    gc()
  }, error = function(e) {
    cat("  ⚠ Hub gene heatmap error:", conditionMessage(e), "\n")
  })
  
  cat("\n  ✅", cell_type_name, "hdWGCNA COMPLETE.\n")
  cat("  Plots saved to:", pdir, "\n\n")
  
  return(seurat_obj)
}

# run hdwgcna for each compartment
# k values are tuned per compartment based on cell counts:
# granulosa: pcos=32,662 + aging=7,890 — k=25 (large, well-represented)
# stromal:   pcos=13,133 + aging=5,108 — k=25 (moderate)
# immune:    pcos=2,316  + aging=2,079 — k=15 (smaller, use lower k)

# load any extra packages needed for plots
if (!requireNamespace("ggrepel",    quietly = TRUE)) install.packages("ggrepel")
if (!requireNamespace("pheatmap",   quietly = TRUE)) install.packages("pheatmap")
if (!requireNamespace("ggpubr",     quietly = TRUE)) install.packages("ggpubr")
library(ggrepel)
library(pheatmap)
library(ggpubr)

# 7a: granulosa cells
# granulosa cells are the primary steroidogenic and folliculogenic cells.
# mitochondrial dysfunction here directly links to impaired folliculogenesis
# in both pcos and aging through reduced atp production and ros accumulation.
cat("\n\n", strrep("█", 70), "\n")
cat("  GRANULOSA CELLS\n")
cat(strrep("█", 70), "\n")

seurat_merged <- run_hdwgcna(
  seurat_obj     = seurat_merged,
  cell_type_name = "Granulosa",
  wgcna_name     = "Granulosa_WGCNA",
  plot_subdir    = "granulosa",
  k_value        = 25,
  fraction       = 0.05,
  min_cells      = 75
)

# save checkpoint after granulosa
saveRDS(seurat_merged, file.path(rds_dir, "seurat_hdwgcna_granulosa.rds"))
cat("✓ Granulosa checkpoint saved.\n")
gc()

# 7b: stromal cells
# stromal cells (fibroblasts, pericytes, theca) regulate the ovarian
# microenvironment. pcos stroma shows androgen excess and fibrosis;
# aged stroma shows senescence-associated secretory phenotype (sasp).
# shared mitochondrial modules here would indicate convergent niche dysfunction.
cat("\n\n", strrep("█", 70), "\n")
cat("  STROMAL CELLS\n")
cat(strrep("█", 70), "\n")

seurat_merged <- run_hdwgcna(
  seurat_obj     = seurat_merged,
  cell_type_name = "Stromal",
  wgcna_name     = "Stromal_WGCNA",
  plot_subdir    = "stromal",
  k_value        = 25,
  fraction       = 0.05,
  min_cells      = 50
)

saveRDS(seurat_merged, file.path(rds_dir, "seurat_hdwgcna_stromal.rds"))
cat("✓ Stromal checkpoint saved.\n")
gc()

# 7c: immune cells
# ovarian immune cells (macrophages, t cells, b cells, nk cells) mediate
# inflammation. both pcos and aging show pro-inflammatory immune activation.
# mitochondrial dysfunction in immune cells drives nlrp3 inflammasome
# activation and impaired oxidative phosphorylation — key shared mechanisms.
cat("\n\n", strrep("█", 70), "\n")
cat("  IMMUNE CELLS\n")
cat(strrep("█", 70), "\n")

seurat_merged <- run_hdwgcna(
  seurat_obj     = seurat_merged,
  cell_type_name = "Immune",
  wgcna_name     = "Immune_WGCNA",
  plot_subdir    = "immune",
  k_value        = 15,       # lower k for smaller immune population
  fraction       = 0.05,
  min_cells      = 30
)

saveRDS(seurat_merged, file.path(rds_dir, "seurat_hdwgcna_all_complete.rds"))
cat("✓ Full object saved.\n")
gc()

# cross-compartment summary
# identify modules shared between pcos and aging across cell types.
# these represent the convergent mitochondrial dysfunction signature.

cat("\n\n", strrep("═", 70), "\n")
cat("  CROSS-COMPARTMENT: Shared Module Summary\n")
cat(strrep("═", 70), "\n\n")

shared_files <- list.files(out_dir, pattern = "shared_modules_", full.names = TRUE)
if (length(shared_files) > 0) {
  shared_all <- do.call(rbind, lapply(shared_files, read.csv))
  cat("Modules significantly dysregulated in BOTH PCOS and Aging:\n")
  print(shared_all)
  write.csv(shared_all,
            file.path(out_dir, "SUMMARY_shared_modules_all_celltypes.csv"),
            row.names = FALSE)
}

# compile enrichment results
enrich_files <- list.files(out_dir, pattern = "enrichment_", full.names = TRUE)
if (length(enrich_files) > 0) {
  enrich_all <- do.call(rbind, lapply(enrich_files, function(f) {
    df <- read.csv(f)
    df$source_file <- basename(f)
    return(df)
  }))
# filter for mito/senescence/ros terms
  mito_terms <- enrich_all %>%
    filter(grepl("mitochond|oxidative|senescen|reactive oxygen|SASP|inflamm|apoptosis|aging",
                 Term, ignore.case = TRUE)) %>%
    arrange(Adjusted.P.value)
  
  write.csv(mito_terms,
            file.path(out_dir, "SUMMARY_mito_senescence_enrichment.csv"),
            row.names = FALSE)
  cat("✓ Mito/senescence enrichment terms extracted:", nrow(mito_terms), "terms.\n")
}

# final summary plot
# one figure summarizing module eigengene activity across all conditions
# and cell types — the "story" figure for the paper.

cat("\nGenerating final summary plot...\n")
tryCatch({
  wgcna_experiments <- c("Granulosa_WGCNA", "Stromal_WGCNA", "Immune_WGCNA")
  cell_type_names   <- c("Granulosa", "Stromal", "Immune")
  
  summary_data <- lapply(seq_along(wgcna_experiments), function(i) {
    wn <- wgcna_experiments[i]
    ct <- cell_type_names[i]
    tryCatch({
      MEs   <- GetMEs(seurat_merged, harmonized = TRUE, wgcna_name = wn)
      meta  <- seurat_merged@meta.data[rownames(MEs), c("condition", "broad_celltype")]
      df    <- cbind(meta, MEs)
      df    <- df %>%
        filter(broad_celltype == ct) %>%
        pivot_longer(cols = -c(condition, broad_celltype),
                     names_to = "Module", values_to = "ME") %>%
        mutate(CellType = ct)
      return(df)
    }, error = function(e) NULL)
  })
  summary_df <- do.call(rbind, Filter(Negate(is.null), summary_data))
  
  if (!is.null(summary_df) && nrow(summary_df) > 0) {
    summary_mean <- summary_df %>%
      group_by(CellType, Module, condition) %>%
      summarise(mean_ME = mean(ME, na.rm = TRUE), .groups = "drop") %>%
      mutate(condition = factor(condition,
                                levels = c("YoungControl","Aged","PCOS_Control","PCOS_Case")))
    
    p_summary <- ggplot(summary_mean,
                        aes(x = condition, y = Module, fill = mean_ME)) +
      geom_tile(color = "white", linewidth = 0.3) +
      facet_wrap(~CellType, scales = "free_y", ncol = 1) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D",
                           midpoint = 0, name = "Mean ME") +
      labs(title    = "Co-expression Module Activity Across Conditions and Cell Types",
           subtitle = "Convergent dysregulation in PCOS and Aging reveals shared mitochondrial programs",
           x = "Condition", y = "Module") +
      theme_cowplot(font_size = 11) +
      theme(axis.text.x  = element_text(angle = 30, hjust = 1),
            strip.text   = element_text(face = "bold", size = 12),
            plot.title    = element_text(face = "bold", hjust = 0.5, size = 13),
            plot.subtitle = element_text(size = 9.5, hjust = 0.5),
            legend.position = "right")
    
    ggsave(file.path(plot_dir, "SUMMARY_module_heatmap_all_conditions.pdf"),
           p_summary, width = 12, height = max(12, length(unique(summary_mean$Module)) * 0.35 + 5),
           dpi = 300)
    ggsave(file.path(plot_dir, "SUMMARY_module_heatmap_all_conditions.png"),
           p_summary, width = 12, height = max(12, length(unique(summary_mean$Module)) * 0.35 + 5),
           dpi = 300)
    cat("✓ Final summary heatmap saved.\n")
  }
}, error = function(e) {
  cat("⚠ Summary plot error:", conditionMessage(e), "\n")
})

# done
cat("\n\n")
cat(strrep("█", 70), "\n")
cat("  ✅ hdWGCNA ANALYSIS COMPLETE\n")
cat(strrep("█", 70), "\n\n")
cat("Output directory:", out_dir, "\n")
cat("Plots directory: ", plot_dir, "\n")
cat("RDS checkpoints: ", rds_dir, "\n\n")
cat("Files generated per compartment (Granulosa / Stromal / Immune):\n")
cat("  01 – Metacell UMAP\n")
cat("  02 – Soft power threshold\n")
cat("  03 – Gene co-expression dendrogram\n")
cat("  04 – Module eigengene feature plots\n")
cat("  05 – Hub gene network\n")
cat("  06 – Module eigengene dot plot by condition\n")
cat("  07 – Module-trait correlation heatmap\n")
cat("  08 – DME volcano (PCOS vs Control, Aged vs Young)\n")
cat("  09 – Enrichment dot plot (GO/KEGG)\n")
cat("  10 – Mito/SASP/Senescence gene set UMAP scores\n")
cat("  11 – Gene set violin plots by condition\n")
cat("  12 – Hub gene expression heatmap\n")
cat("  SUMMARY – Cross-compartment module heatmap\n\n")

# session info for reproducibility
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo.txt"))
cat("Session info saved to:", file.path(out_dir, "sessionInfo.txt"), "\n")