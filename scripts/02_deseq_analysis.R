source("scripts/setup.R")

# Loading data
load("results/QC/Filtered_Data.RData")

# -------------------------------------
# DESEQ2 DIFFERENTIAL EXPRESSION ANALYSIS
# -------------------------------------

cat("\nDIFFERENTIAL EXPRESSION ANALYSIS - DESeq2\n")

dds <- DESeqDataSetFromMatrix(countData = count_matrix_filtered,
                              colData = sample_info,
                              design = ~ Condition)
cat("Running DESeq2 analysis...\n")
dds <- DESeq(dds)

cat("\nAvailable contrasts:\n")
print(resultsNames(dds))
cat("\n")

# COMPARISON 1: Resistant vs Sensitive (Control)
cat("=== COMPARISON 1: Resistant vs Sensitive (Control) ===\n")

res_resistant_vs_sensitive <- results(dds, 
                                      contrast = c("Condition",
                                                   "Resistant_Control",
                                                   "Sensitive_Control"),
                                      alpha = 0.05)

summary(res_resistant_vs_sensitive, alpha = 0.05)

res_resistant_vs_sensitive_ordered <- res_resistant_vs_sensitive[order(res_resistant_vs_sensitive$padj), ]

sig_genes_comp1 <- sum(res_resistant_vs_sensitive$padj < 0.05, na.rm = TRUE)
up_genes_comp1 <- sum(res_resistant_vs_sensitive$padj < 0.05 & 
                        res_resistant_vs_sensitive$log2FoldChange > 1, na.rm = TRUE)
down_genes_comp1 <- sum(res_resistant_vs_sensitive$padj < 0.05 & 
                          res_resistant_vs_sensitive$log2FoldChange < -1, na.rm = TRUE)

cat("\nStatistics:\n")
cat("  Significant genes (padj < 0.05):", sig_genes_comp1, "\n")
cat("  Upregulated (|log2FC| > 1):", up_genes_comp1, "\n")
cat("  Downregulated (|log2FC| > 1):", down_genes_comp1, "\n\n")

cat("Top 10 DEGs:\n")
print(head(res_resistant_vs_sensitive_ordered[, c("baseMean", "log2FoldChange", "padj")], 10))

write.csv(as.data.frame(res_resistant_vs_sensitive_ordered),
          "results/DESeq2/DESeq2_Resistant_vs_Sensitive_Control.csv", row.names = TRUE)
cat("\nSaved: results/DESeq2/DESeq2_Resistant_vs_Sensitive_Control.csv\n\n")

# --- COMPARISON 2: Brequinar Effect in Resistant Cells ---
cat("=== COMPARISON 2: Brequinar in Resistant Cells ===\n")

res_resistant_BQ_vs_Control <- results(dds, 
                                       contrast = c("Condition",
                                                    "Resistant_Brequinar",
                                                    "Resistant_Control"),
                                       alpha = 0.05)

summary(res_resistant_BQ_vs_Control, alpha = 0.05)

res_resistant_BQ_vs_Control_ordered <- res_resistant_BQ_vs_Control[order(res_resistant_BQ_vs_Control$padj), ]

sig_genes_comp2 <- sum(res_resistant_BQ_vs_Control$padj < 0.05, na.rm = TRUE)
up_genes_comp2 <- sum(res_resistant_BQ_vs_Control$padj < 0.05 & 
                        res_resistant_BQ_vs_Control$log2FoldChange > 1, na.rm = TRUE)
down_genes_comp2 <- sum(res_resistant_BQ_vs_Control$padj < 0.05 & 
                          res_resistant_BQ_vs_Control$log2FoldChange < -1, na.rm = TRUE)

cat("\nStatistics:\n")
cat("  Significant genes (padj < 0.05):", sig_genes_comp2, "\n")
cat("  Upregulated (|log2FC| > 1):", up_genes_comp2, "\n")
cat("  Downregulated (|log2FC| > 1):", down_genes_comp2, "\n\n")

cat("Top 10 DEGs:\n")
print(head(res_resistant_BQ_vs_Control_ordered[, c("baseMean", "log2FoldChange", "padj")], 10))

write.csv(as.data.frame(res_resistant_BQ_vs_Control_ordered),
          "results/DESeq2/DESeq2_Resistant_Brequinar_vs_Control.csv", row.names = TRUE)
cat("\nSaved: results/DESeq2/DESeq2_Resistant_Brequinar_vs_Control.csv\n\n")

# --- COMPARISON 3: Brequinar Effect in Sensitive Cells ---
cat("=== COMPARISON 3: Brequinar in Sensitive Cells ===\n")

res_sensitive_BQ_vs_Control <- results(dds, 
                                       contrast = c("Condition",
                                                    "Sensitive_Brequinar",
                                                    "Sensitive_Control"),
                                       alpha = 0.05)

summary(res_sensitive_BQ_vs_Control, alpha = 0.05)

res_sensitive_BQ_vs_Control_ordered <- res_sensitive_BQ_vs_Control[order(res_sensitive_BQ_vs_Control$padj), ]

sig_genes_comp3 <- sum(res_sensitive_BQ_vs_Control$padj < 0.05, na.rm = TRUE)
up_genes_comp3 <- sum(res_sensitive_BQ_vs_Control$padj < 0.05 & 
                        res_sensitive_BQ_vs_Control$log2FoldChange > 1, na.rm = TRUE)
down_genes_comp3 <- sum(res_sensitive_BQ_vs_Control$padj < 0.05 & 
                          res_sensitive_BQ_vs_Control$log2FoldChange < -1, na.rm = TRUE)

cat("\nStatistics:\n")
cat("  Significant genes (padj < 0.05):", sig_genes_comp3, "\n")
cat("  Upregulated (|log2FC| > 1):", up_genes_comp3, "\n")
cat("  Downregulated (|log2FC| > 1):", down_genes_comp3, "\n\n")

cat("Top 10 DEGs:\n")
print(head(res_sensitive_BQ_vs_Control_ordered[, c("baseMean", "log2FoldChange", "padj")], 10))

write.csv(as.data.frame(res_sensitive_BQ_vs_Control_ordered),
          "results/DESeq2/DESeq2_Sensitive_Brequinar_vs_Control.csv", row.names = TRUE)
cat("\nSaved: results/DESeq2/DESeq2_Sensitive_Brequinar_vs_Control.csv\n\n")

# --- COMPARISON 4: Interaction Effect ---
cat("=== COMPARISON 4: Interaction Effect ===\n")

dds_interaction <- DESeqDataSetFromMatrix(countData = count_matrix_filtered,
                                          colData = sample_info,
                                          design = ~ Resistance + Treatment + Resistance:Treatment)

cat("Running interaction model...\n")
dds_interaction <- DESeq(dds_interaction)

cat("\nAvailable interaction results:\n")
print(resultsNames(dds_interaction))


res_interaction <- results(dds_interaction, 
                           name = "ResistanceSensitive.TreatmentControl",  
                           alpha = 0.05)

summary(res_interaction, alpha = 0.05)

res_interaction_ordered <- res_interaction[order(res_interaction$padj), ]

sig_genes_comp4 <- sum(res_interaction$padj < 0.05, na.rm = TRUE)
up_genes_comp4 <- sum(res_interaction$padj < 0.05 & 
                        res_interaction$log2FoldChange > 1, na.rm = TRUE)
down_genes_comp4 <- sum(res_interaction$padj < 0.05 & 
                          res_interaction$log2FoldChange < -1, na.rm = TRUE)

cat("\nStatistics:\n")
cat("Significant interaction genes:", sig_genes_comp4, "\n")
cat("Positive interaction (|log2FC| > 1):", up_genes_comp4, "\n")
cat("Negative interaction (|log2FC| > 1):", down_genes_comp4, "\n\n")

cat("Interpretation:\n")
cat("Positive log2FC: Effect is STRONGER in Sensitive+Control combination\n")
cat("Negative log2FC: Effect is WEAKER in Sensitive+Control combination\n\n")

if(sig_genes_comp4 > 0) {
  cat("Top 10 Genes with Differential Response:\n")
  print(head(res_interaction_ordered[, c("baseMean", "log2FoldChange", "padj")], 10))
} else {
  cat("No significant interaction effects found\n")
}

write.csv(as.data.frame(res_interaction_ordered),
          "results/DESeq2/DESeq2_Interaction_Effect.csv", row.names = TRUE)
cat("\nSaved: results/DESeq2/DESeq2_Interaction_Effect.csv\n\n")

# --- SUMMARY TABLE ---

cat("SUMMARY OF ALL COMPARISONS\n")

summary_stats <- data.frame(
  Comparison = c("Resistant vs Sensitive (Control)",
                 "Brequinar in Resistant Cells",
                 "Brequinar in Sensitive Cells",
                 "Interaction (Differential Response)"),
  Total_Genes = rep(nrow(count_matrix_filtered), 4),
  Sig_DEGs = c(sig_genes_comp1, sig_genes_comp2, sig_genes_comp3, sig_genes_comp4),
  Upregulated = c(up_genes_comp1, up_genes_comp2, up_genes_comp3, up_genes_comp4),
  Downregulated = c(down_genes_comp1, down_genes_comp2, down_genes_comp3, down_genes_comp4),
  Percent_Sig = round(c(sig_genes_comp1, sig_genes_comp2, sig_genes_comp3, sig_genes_comp4)/
                        nrow(count_matrix_filtered)*100, 2)
)

print(summary_stats)

write.csv(summary_stats, "results/DESeq2/DESeq2_Analysis_Summary.csv", row.names = FALSE)
cat("\nSaved: results/DESeq2/DESeq2_Analysis_Summary.csv\n\n")

# --- SAVE R OBJECTS ---
save(dds, dds_interaction, 
     res_resistant_vs_sensitive_ordered,
     res_resistant_BQ_vs_Control_ordered,
     res_sensitive_BQ_vs_Control_ordered,
     res_interaction_ordered,
     sample_info,
     summary_stats,
     file = "results/DESeq2/DESeq2_Results.RData")

cat("All DESeq2 objects saved: results/DESeq2/DESeq2_Results.RData\n")