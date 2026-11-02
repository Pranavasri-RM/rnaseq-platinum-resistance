# Packages used across all steps

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tibble)
  library(ggrepel)
})

# Make output folders so the scripts also run outside Snakemake
for (d in c("results/QC", "results/DESeq2", "results/visualization",
            "results/enrichment", "results/biomarkers", "results/figures")) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}