#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# usage: rscript stage4_sct_integrate.r <dataset>
# dataset ? {aging, pcos}

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("usage: Rscript stage4_sct_integrate.R <dataset>")
dataset <- args[1]

# directories
base_output <- "~/projects/scrna/output"

in_dir       <- file.path(base_output, "03_doublets", dataset)
per_sample_out <- file.path(base_output, "04_normalization", dataset)
integrate_out  <- file.path(base_output, "05_integration", dataset)
plots_out      <- file.path(integrate_out, "plots")

dir.create(per_sample_out, recursive = TRUE, showWarnings = FALSE)
dir.create(integrate_out,  recursive = TRUE, showWarnings = FALSE)
dir.create(plots_out,      recursive = TRUE, showWarnings = FALSE)

# integration parameters
sct_vars_to_regress   <- "percent.mt"
select_nfeatures      <- 3000
integration_dims       <- 1:30
integration_resolution <- 0.5

# detect sample files
pattern <- paste0("^", dataset, "_sample_.*\\.rds$")
sample_files <- list.files(in_dir, pattern = pattern, full.names = TRUE)

if (length(sample_files) == 0)
  stop("no standardized sample files found in: ", in_dir)

message("Found ", length(sample_files), " sample files. Running SCTransform per sample...\n")

# per-sample normalization (fresh sctransform models)
slist <- list()

for (f in sample_files) {

# extract sample suffix safely
  samp <- gsub("\\.rds$", "", basename(f))
  samp <- sub(paste0(dataset, "_sample_"), "", samp)

  message(" loading: ", basename(f))
  so <- readRDS(f)
  DefaultAssay(so) <- "RNA"

# collapse rna layers if needed
  if (length(so[["RNA"]]@layers) > 1) {
    message("  joining RNA layers")
    so <- JoinLayers(so, assay = "RNA")
  }

# ensure percent.mt exists
  if (is.null(so$percent.mt) || all(is.na(so$percent.mt))) {
    so$percent.mt <- PercentageFeatureSet(so, pattern="^mt|^MT|^mt-")
  }

  message("  running SCTransform on: ", samp)
  so <- SCTransform(so,
                    vars.to.regress = sct_vars_to_regress,
                    verbose = FALSE)

  outpath <- file.path(per_sample_out, paste0(samp, "_sct.rds"))
  saveRDS(so, outpath)
  slist[[samp]] <- so

  message("  saved: ", outpath, "\n")
}

# integration
message("Selecting integration features")
features <- SelectIntegrationFeatures(slist, nfeatures = select_nfeatures)

message("Prep SCT Integration")
slist <- PrepSCTIntegration(object.list = slist,
                            anchor.features = features,
                            verbose = TRUE)

message("Finding anchors")
anchors <- FindIntegrationAnchors(object.list = slist,
                                  normalization.method = "SCT",
                                  anchor.features = features,
                                  verbose = TRUE)

message("Integrating data")
combined <- IntegrateData(anchorset = anchors,
                          normalization.method = "SCT",
                          verbose = TRUE)

# dimensionality reduction
message("running PCA")
combined <- RunPCA(combined, verbose = FALSE)

message("running UMAP")
combined <- RunUMAP(combined, dims = integration_dims, verbose = FALSE)

combined <- FindNeighbors(combined, dims = integration_dims, verbose = FALSE)
combined <- FindClusters(combined,
                         resolution = integration_resolution,
                         verbose = FALSE)

# save integrated object
integrated_file <- file.path(integrate_out,
                             paste0(dataset, "_integrated.rds"))
saveRDS(combined, integrated_file)
message("saved integrated object: ", integrated_file)

# umap plot
p <- DimPlot(combined, reduction = "umap", label = TRUE) +
     ggtitle(paste0(dataset, " integrated UMAP"))

ggsave(file.path(plots_out, "umap_clusters.pdf"),
       p, width = 8, height = 6)

# save integration parameters
params <- data.frame(
  name = c("sct_vars_to_regress",
           "select_nfeatures",
           "integration_dims",
           "integration_resolution"),
  value = c(sct_vars_to_regress,
            select_nfeatures,
            paste(integration_dims, collapse=","),
            integration_resolution),
  stringsAsFactors = FALSE
)

write.csv(params,
          file.path(integrate_out, "integration_parameters.csv"),
          row.names = FALSE)

message("\nstage4 complete. outputs -> ", integrate_out, "\n")