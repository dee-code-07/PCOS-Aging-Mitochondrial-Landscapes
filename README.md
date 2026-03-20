# Single-Cell and Spatial Landscapes of Mitochondrial Dysfunction in PCOS and Ovarian Aging

## Overview
This repository contains the official R implementation for the multi-modal analysis of mitochondrial dysfunction in Polycystic Ovary Syndrome (PCOS) and Ovarian Aging. We employ integrated single-cell RNA sequencing (scRNA-seq) and Spatial Transcriptomics to map the cellular and spatial dynamics of mitochondrial decline across these two conditions.

---

## Repository Structure

The scripts are organized sequentially to reflect the workflow presented in the manuscript:

### 📂 `01_preprocessing/`
Clean, step-by-step pipelines for initial data processing.
- **`scrna/`**: Quality control, doublet removal (scDblFinder), normalization, Harmony integration, and cluster annotation.
- **`spatial/`**: Spatial-specific QC, SCTransform normalization with custom PCOS-fix, and label transfer from scRNA-seq references.

### 📂 `02_analysis/`
The core analysis scripts used to generate primary figures:
- **`fig1_overview`**: UMAP visualizations and cell-type composition analysis.
- **`fig2_modules`**: Mitochondrial and inflammatory module scoring and statistical modeling.
- **`fig3_deg_enrichment`**: Differential expression and functional enrichment (GO/KEGG).
- **`fig4_hdwgcna`**: High-dimensional Weighted Gene Co-expression Network Analysis.
- **`fig5_pseudotime`**: Trajectory inference and lineage analysis.
- **`fig6_spatial`**: Spatially variable gene (SVG) detection and spatial cell-type mapping.

### 📂 `03_supplementary/`
Additional validation and supplementary analyses including Venn overlap, gene prioritization, network analysis, and transcription factor (TF) enrichment.

---

## Dependencies
To run these scripts, you will need **R (>= 4.2.0)** and the following key packages:
- **Core**: `Seurat (v4/v5)`, `tidyverse`, `patchwork`, `ggplot2`
- **Network**: `hdWGCNA`, `WGCNA`, `igraph`
- **Trajectory**: `monocle3`, `SeuratWrappers`
- **Statistics**: `lme4`, `ggpubr`, `clusterProfiler`
- **Spatial**: `Giotto` or `Seurat` spatial modules

---

## Usage
Scripts are numbered sequentially. For best results, follow the order:
1. Run `01_preprocessing` to generate integrated objects.
2. Execute `02_analysis` folders to reproduce manuscript figures.
3. Consult `03_supplementary` for deeper metric validation.

---

## Citation
If you use these scripts or find this research helpful, please cite:
> *Manuscript Title*, Authors, Journal (Year). DOI: [Link]

## Contact
For questions regarding the code or data, please open an issue or contact the authors directly.
