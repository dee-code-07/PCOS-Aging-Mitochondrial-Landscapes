#!/usr/bin/env Rscript

# disable parallel workers for findintegrationanchors / integratedata
suppressPackageStartupMessages(library(future))
plan("sequential")
options(future.globals.maxSize = Inf)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

# user config
sct_base <- "~/projects/spatial/output/03_sct"   # stage2 output
out_base <- "~/projects/spatial/output/04_integration"
dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

datasets <- c("pcos", "aging")   # names of dataset folders under sct_base
n_dims <- 30
cluster_resolution <- 0.4

# whether to run joint integration across datasets (pcos+aging)
do_joint_integration <- FALSE

# helper to write plots
safe_save_plot <- function(plot_obj, file, width=7, height=6) {
  tryCatch({
    ggsave(file, plot_obj, width=width, height=height, dpi=300)
  }, error = function(e) message("Could not save plot: ", file, " : ", e$message))
}

# functions
load_sct_samples <- function(dataset_name) {
  dir_samples <- file.path(sct_base, dataset_name)
  if (!dir.exists(dir_samples)) {
    message("No SCT directory for ", dataset_name, " at ", dir_samples)
    return(list())
  }
  rds_files <- list.files(dir_samples, pattern = "_sct\\.rds$", full.names = TRUE)
  objs <- lapply(rds_files, function(f) {
    message("  loading ", f)
    readRDS(f)
  })
  names(objs) <- tools::file_path_sans_ext(basename(rds_files))
  return(objs)
}

run_integration_on_list <- function(obj_list, out_prefix, outdir) {
# obj_list: named list of seurat objects already sctransformed (sct assay present)
  if (length(obj_list) == 0) {
    message("No objects supplied for ", out_prefix)
    return(NULL)
  }

# ensure sct is active and used
  for (nm in names(obj_list)) {
    so <- obj_list[[nm]]
    DefaultAssay(so) <- "SCT"
    obj_list[[nm]] <- so
  }

  message("Preparing features for integration...")
  features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)
  obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features, verbose = FALSE)

  message("Finding anchors (SCT)...")
  anchors <- FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT",
                                    anchor.features = features, verbose = FALSE)

  message("Integrating data for ", out_prefix)
  int <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)

  DefaultAssay(int) <- "integrated"  # integrated assay
# run standard workflow on integrated assay
  int <- ScaleData(int, verbose = FALSE)
  int <- RunPCA(int, verbose = FALSE, npcs = n_dims)
  int <- RunUMAP(int, reduction = "pca", dims = 1:n_dims, verbose = FALSE)
  int <- FindNeighbors(int, reduction = "pca", dims = 1:n_dims, verbose = FALSE)
  int <- FindClusters(int, resolution = cluster_resolution, verbose = FALSE)

# save integrated object + plots
  out_rds <- file.path(outdir, paste0(out_prefix, "_integrated.rds"))
  saveRDS(int, out_rds)
  message("Saved integrated object: ", out_rds)

# plots
  p1 <- DimPlot(int, reduction = "umap", group.by = "orig.ident") + ggtitle(paste0(out_prefix, " UMAP by sample"))
  p2 <- DimPlot(int, reduction = "umap", label = TRUE) + ggtitle(paste0(out_prefix, " UMAP clusters"))
  safe_save_plot(p1, file.path(outdir, paste0(out_prefix, "_UMAP_by_sample.png")))
  safe_save_plot(p2, file.path(outdir, paste0(out_prefix, "_UMAP_clusters.png")))

# cluster sizes
  cluster_counts <- table(Idents(int))
  write.csv(as.data.frame(cluster_counts), file.path(outdir, paste0(out_prefix, "_cluster_counts.csv")), row.names = FALSE)

# find markers per cluster (top 5)
  message("Finding markers for clusters (integrated)...")
  DefaultAssay(int) <- "integrated"
  markers <- tryCatch({
    FindAllMarkers(int, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
  }, error = function(e) {
    message("FindAllMarkers failed: ", e$message)
    NULL
  })
  if (!is.null(markers)) {
    write.csv(markers, file.path(outdir, paste0(out_prefix, "_markers_allclusters.csv")), row.names = FALSE)
# top markers per cluster
    top_markers <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 5)
    write.csv(top_markers, file.path(outdir, paste0(out_prefix, "_top5_markers.csv")), row.names = FALSE)
  }

  return(int)
}

# main
message("Stage 3 — Integration starting")

dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

# per-dataset integration
integrated_objects <- list()
for (ds in datasets) {
  message("\nProcessing dataset: ", ds)
  objs <- load_sct_samples(ds)
  if (length(objs) <= 1) {
    message("  Only one or zero samples for ", ds, " — skipping integration; saving single object copy instead.")
    if (length(objs) == 1) {
      saveRDS(objs[[1]], file.path(out_base, paste0(ds, "_single_sct_copy.rds")))
      integrated_objects[[ds]] <- objs[[1]]
    }
    next
  }
  outdir <- file.path(out_base, ds)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  int_obj <- run_integration_on_list(objs, out_prefix = ds, outdir = outdir)
  integrated_objects[[ds]] <- int_obj
}

# optional joint integration across datasets (pcos + aging)
if (do_joint_integration) {
  message("\nStarting joint integration (pcos + aging). This may be memory intensive.")
# use integrated per-dataset objects if present, else use sample lists directly
  joint_list <- list()
  if (!is.null(integrated_objects$pcos) && inherits(integrated_objects$pcos, "Seurat")) {
# split integrated back into samples by orig.ident
    pcos_int <- integrated_objects$pcos
    for (s in unique(pcos_int$orig.ident)) {
      sub <- subset(pcos_int, subset = orig.ident == s)
      DefaultAssay(sub) <- "SCT"  # ensure SCT if available; otherwise we will re-run SCT later if needed
      joint_list[[paste0("pcos_", s)]] <- sub
    }
  } else {
# fallback: load raw sct samples
    objs_p <- load_sct_samples("pcos")
    joint_list <- c(joint_list, objs_p)
  }

  if (!is.null(integrated_objects$aging) && inherits(integrated_objects$aging, "Seurat")) {
    aging_int <- integrated_objects$aging
    for (s in unique(aging_int$orig.ident)) {
      sub <- subset(aging_int, subset = orig.ident == s)
      DefaultAssay(sub) <- "SCT"
      joint_list[[paste0("aging_", s)]] <- sub
    }
  } else {
    objs_a <- load_sct_samples("aging")
    joint_list <- c(joint_list, objs_a)
  }

# finally run integration on joint_list
  if (length(joint_list) >= 2) {
    joint_outdir <- file.path(out_base, "joint_pcos_aging")
    dir.create(joint_outdir, recursive = TRUE, showWarnings = FALSE)
    joint_int <- run_integration_on_list(joint_list, out_prefix = "pcos_aging_joint", outdir = joint_outdir)
    integrated_objects[["joint_pcos_aging"]] <- joint_int
  } else {
    message("Not enough objects to run joint integration.")
  }
} # end joint integration

# save a small summary
summary_rows <- list()
for (nm in names(integrated_objects)) {
  obj <- integrated_objects[[nm]]
  if (is.null(obj)) next
  ncells <- ncol(obj)
  nclusters <- length(unique(Idents(obj)))
  summary_rows[[length(summary_rows) + 1]] <- data.frame(name = nm, cells = ncells, nclusters = nclusters)
}
if (length(summary_rows) > 0) {
  sumdf <- do.call(rbind, summary_rows)
  write.csv(sumdf, file.path(out_base, "integration_summary_stage3.csv"), row.names = FALSE)
  message("Integration summary saved.")
}

message("Stage 3 complete.")
