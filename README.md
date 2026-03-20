# Single-cell and spatial landscape of mitochondrial dysfunction in PCOS and ovarian aging

**Authors:** Deeksha H and Budheswar Dehury*  
**Affiliation:** Department of Bioinformatics, Manipal School of Life Sciences, Manipal Academy of Higher Education, Manipal-576104, India  
**Correspondence:** budheswar.dehury@manipal.edu

---

## Overview
This repository contains the official R implementation for the multi-modal analysis of mitochondrial dysfunction in Polycystic Ovary Syndrome (PCOS) and Ovarian Aging. By integrating mouse single-cell RNA-seq and Visium spatial transcriptomics, this study maps the cellular and spatial dynamics of mitochondrial decline across these two conditions.

### Key Findings
- **Accelerated Molecular Aging:** Granulosa cells in PCOS occupy a transcriptional state intermediate between young and aged cells, suggesting accelerated ovarian aging at single-cell resolution.
- **Convergent Markers:** Identified 83 concordantly dysregulated shared genes and 17 high-confidence priority candidates (*Gstp1*, *Gja1*, *Fst*).
- **Mitochondrial-Senescence Disconnect:** PCOS presents a unique disconnect where high bioenergetic stress coincides with a blunted physiological senescence program required for luteinization.
- **Divergent OXPHOS Trajectories:** High-dimensional co-expression network analysis (hdWGCNA) uncovered that Oxidative Phosphorylation (OXPHOS) is suppressed in PCOS but upregulated in aging, pointing to distinct routes to follicular failure.

---

## Data Availability
The analysis uses publicly available datasets from the Gene Expression Omnibus (GEO):
- **PCOS scRNA-seq:** [GSE268919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268919)
- **Ovarian Aging scRNA-seq:** [GSE232309](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE232309)
- **PCOS Spatial Transcriptomics:** [GSE296728](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE296728)
- **Ovarian Aging Spatial Transcriptomics:** [GSE188257](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188257)

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
