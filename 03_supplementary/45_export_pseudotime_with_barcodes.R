suppressPackageStartupMessages(library(Seurat))
obj <- readRDS("/home/deekshah/mini_project/analysis/29_joint_pseudotime_rpca/merged_rpca_harmony_checkpoint.rds")
meta <- obj@meta.data
meta$cell_barcode <- rownames(meta)
write.csv(meta, "/home/deekshah/mini_project/revision_final/pseudotime_data_with_barcodes.csv", row.names=FALSE)
