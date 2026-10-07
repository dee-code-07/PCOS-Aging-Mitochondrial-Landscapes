# modified script with unified y-axis
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

OUTPUT_DIR <- "E:/Documents/mini_project/analysis/04_module_statistics"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

pcos_meta <- read.csv("/mnt/e/mini_project/revision_final/meta_task12.csv") %>% filter(dataset == "PCOS")

modules <- list(
  Mitochondrial = "Mitochondrial_Score1",
  Senescence    = "Senescence_Score1",
  Oxidative     = "Oxidative_Score1",
  Inflammation  = "Inflammation_Score1"
)

all_scores <- unlist(pcos_meta[, unname(unlist(modules))])
g_min <- min(all_scores, na.rm=TRUE)
g_max <- max(all_scores, na.rm=TRUE)
# Add some padding for the p-value brackets
g_max_pad <- g_max + (g_max - g_min) * 0.15

cat(sprintf("Exact y-axis range used: [%f, %f]\n", g_min, g_max_pad))

theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.text    = element_text(colour = "black"),
      axis.title.x = element_blank(),
      plot.title   = element_text(face = "bold", hjust = 0.5, size = base_size + 1),
      legend.position = "none"
    )
}

pcos_fills  <- c(case = "#D55E00", control = "#0072B2")
# Wait, what are the conditions in the saved metadata? 
# In task12_metadata.R, it's just meta_p$condition!
