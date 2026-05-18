source("scripts/setup.R")

# From DESeq step
load("results/DESeq2/DESeq2_Results.RData")

# From annotation step
load("results/DESeq2/DESeq2_Annotated_Results.RData")

# From pyrimidine step
load("results/enrichment/Pyrimidine_Metabolism_Results.RData")

# ============================================
# IDENTIFY BIOMARKER CANDIDATES
# ============================================

cat("\nBIOMARKER CANDIDATE IDENTIFICATION\n")

# Biomarker selection criteria
biomarker_candidates <- res_resistant_annotated %>%
  filter(padj < 0.01,
         abs(log2FoldChange) > 2,
         baseMean > 100) %>%
  arrange(padj)

cat("Found", nrow(biomarker_candidates), "biomarker candidates\n\n")

# Getting normalized counts
biomarker_genes <- biomarker_candidates$gene_symbol[1:20]
biomarker_counts <- counts(dds, normalized = TRUE)[biomarker_genes, ]

# Plotting top biomarkers
plot_biomarker <- function(gene_name, counts_matrix, sample_info) {
  gene_expr <- as.numeric(counts_matrix[gene_name, ])
  
  plot_df <- data.frame(
    Expression = log2(gene_expr + 1),
    Resistance = sample_info$Resistance,
    Treatment = sample_info$Treatment
  )
  
  ggplot(plot_df, aes(x = Resistance, y = Expression, fill = Treatment)) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge()) +
    theme_bw() +
    labs(title = gene_name, y = "log2(Normalized Counts + 1)") +
    scale_fill_manual(values = c("Control" = "lightblue", "Brequinar" = "coral"))
}

pdf("results/biomarkers/Biomarker_Candidates_Expression.pdf", width = 12, height = 10)
for (i in 1:min(12, nrow(biomarker_candidates))) {
  gene <- biomarker_candidates$gene_symbol[i]
  if (!is.na(gene) && gene %in% rownames(biomarker_counts)) {
    print(plot_biomarker(gene, biomarker_counts, sample_info))
  }
}
dev.off()

cat("Saved: results/biomarkers/Biomarker_Candidates_Expression.pdf\n\n")

# Save results
write.csv(biomarker_candidates,
          "results/biomarkers/Biomarker_Candidates_for_Brequinar_Resistance.csv",
          row.names = FALSE)

cat("Saved: results/biomarkers/Biomarker_Candidates_for_Brequinar_Resistance.csv\n")

# Summary
cat("\nBiomarker Summary:\n")
cat("Total candidates:", nrow(biomarker_candidates), "\n")
cat("Upregulated:", sum(biomarker_candidates$log2FoldChange > 2), "\n")
cat("Downregulated:", sum(biomarker_candidates$log2FoldChange < -2), "\n\n")

save(biomarker_candidates, file = "results/biomarkers/Biomarker_Results.RData")
cat("All results saved: results/biomarkers/Biomarker_Results.RData\n")

# ============================================
# FINAL PIPELINE SUMMARY
# ============================================

cat("FINAL PIPELINE SUMMARY\n")

# Basic counts

total_genes <- nrow(res_resistant_annotated)

sig_genes <- sum(res_resistant_annotated$padj < 0.05, na.rm = TRUE)

biomarker_count <- nrow(biomarker_candidates)

pyrimidine_total <- sum(
  ifelse(is.null(pyrimidine_resistant), 0, nrow(pyrimidine_resistant)),
  ifelse(is.null(pyrimidine_BQ_resistant), 0, nrow(pyrimidine_BQ_resistant)),
  ifelse(is.null(pyrimidine_BQ_sensitive), 0, nrow(pyrimidine_BQ_sensitive)),
  ifelse(is.null(pyrimidine_interaction), 0, nrow(pyrimidine_interaction))
)

# Creating summary table

final_summary <- data.frame(
  Metric = c(
    "Total Genes Analyzed",
    "Significant DEGs (padj < 0.05)",
    "Biomarker Candidates (padj < 0.01 & |log2FC| > 2)",
    "Total Pyrimidine-related Genes Identified"
  ),
  Value = c(
    total_genes,
    sig_genes,
    biomarker_count,
    pyrimidine_total
  )
)

print(final_summary)

# Saving summary

write.csv(final_summary,
          "results/biomarkers/Final_Pipeline_Summary.csv",
          row.names = FALSE)

cat("Saved: results/biomarkers/Final_Pipeline_Summary.csv\n")
