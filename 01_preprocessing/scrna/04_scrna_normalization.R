#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tools)   # for file_path_sans_ext
})

# usage: rscript stage3_no_sct.r <dataset> <rename_mode>
# rename_mode: "b" to rename to standardized names (recommended), "a" to keep original
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("usage: Rscript stage3_no_sct.R <dataset> [rename_mode]")
dataset <- args[1]
rename_mode <- ifelse(length(args)>=2, args[2], "b")

base_output <- "~/projects/scrna/output"
in_dir <- file.path(base_output, "03_doublets", dataset)
out_dir <- in_dir  # we keep files in same folder but may rename
summary_csv <- file.path(out_dir, "stage3_summary.csv")

files <- list.files(in_dir, pattern = "_doubletfiltered\\.rds$", full.names = TRUE)
if (length(files)==0) stop("no doubletfiltered files found in: ", in_dir)

summary <- data.frame(sample_orig=character(), sample_new=character(), path_orig=character(), path_new=character(),
                      cells=integer(), stringsAsFactors=FALSE)

i <- 1
for (f in files) {
  fname <- basename(f)
  samp_orig <- sub("_doubletfiltered\\.rds$", "", fname)
  so <- tryCatch(readRDS(f), error=function(e) { warning("cannot read: ", f); NULL })
  if (is.null(so)) next
  cells <- ncol(so)

  if (tolower(rename_mode) == "b") {
# standardized name: <dataset>_sample_<basename>.rds
    samp_new <- paste0(dataset, "_sample_", samp_orig)
    newpath <- file.path(out_dir, paste0(samp_new, ".rds"))
# save a copy with new name (do not remove original)
    saveRDS(so, newpath)
  } else {
    samp_new <- samp_orig
    newpath <- f
  }

  summary[i,] <- list(samp_orig, samp_new, f, newpath, cells)
  i <- i + 1
  message("processed: ", samp_orig, " -> ", samp_new, " (n=", cells, ")")
}

write.csv(summary, summary_csv, row.names = FALSE)
message("stage3 complete. summary -> ", summary_csv)
