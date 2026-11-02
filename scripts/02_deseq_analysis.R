source("scripts/setup.R")

load("results/QC/Filtered_Data.RData")

# -------------------------------------
# DIFFERENTIAL EXPRESSION (DESeq2)
#   Two models on the same counts:
#   - ~ Condition for the three comparisons: resistant vs sensitive (untreated),
#     and brequinar vs control within each line
#   - ~ Resistance * Treatment to test whether brequinar acts differently in
#     resistant cells (the interaction term)
#   With 3 replicates, fold changes for low-count genes are noisy, so I shrink
#   them with ashr. This changes log2FoldChange and lfcSE only, not p-values.
# -------------------------------------

cat("\nDIFFERENTIAL EXPRESSION ANALYSIS - DESeq2\n")

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1

dds <- DESeqDataSetFromMatrix(countData = count_matrix_filtered,
                              colData   = sample_info,
                              design    = ~ Condition)
cat("Running DESeq2 (main model)...\n")
dds <- DESeq(dds)

cat("\nCoefficients:\n")
print(resultsNames(dds))
cat("\n")

# -------------------------------------
# Helper: prints the summary and top genes for one comparison, saves the CSV.
#   I report counts with and without |log2FC| > 1, since both get quoted.
# -------------------------------------

summarise_and_save <- function(res, label, csv_path) {
  cat("===", label, "===\n")
  summary(res, alpha = PADJ_CUTOFF)

  res_ordered <- res[order(res$padj, na.last = TRUE), ]

  sig  <- !is.na(res$padj) & res$padj < PADJ_CUTOFF
  stats <- c(
    Sig_padj05     = sum(sig),
    Up_padj05      = sum(sig & res$log2FoldChange > 0),
    Down_padj05    = sum(sig & res$log2FoldChange < 0),
    Up_padj05_LFC1   = sum(sig & res$log2FoldChange >  LFC_CUTOFF),
    Down_padj05_LFC1 = sum(sig & res$log2FoldChange < -LFC_CUTOFF)
  )

  cat("\n  padj < 0.05:                ", stats["Sig_padj05"],
      " (", stats["Up_padj05"], " up, ", stats["Down_padj05"], " down)\n", sep = "")
  cat("  padj < 0.05 and |log2FC| > 1:", stats["Up_padj05_LFC1"], "up,",
      stats["Down_padj05_LFC1"], "down\n\n")

  cat("Top 10 genes by padj:\n")
  print(head(res_ordered[, c("baseMean", "log2FoldChange", "padj")], 10))

  write.csv(as.data.frame(res_ordered), csv_path, row.names = TRUE)
  cat("\nSaved:", csv_path, "\n\n")

  list(res = res_ordered, stats = stats)
}

# -------------------------------------
# COMPARISONS 1-3 (main model)
# -------------------------------------

shrink_contrast <- function(contrast) {
  res <- results(dds, contrast = contrast, alpha = PADJ_CUTOFF)
  lfcShrink(dds, contrast = contrast, res = res, type = "ashr", quiet = TRUE)
}

comp1 <- summarise_and_save(
  shrink_contrast(c("Condition", "Resistant_Control", "Sensitive_Control")),
  "COMPARISON 1: Resistant vs Sensitive (untreated)",
  "results/DESeq2/DESeq2_Resistant_vs_Sensitive_Control.csv")

comp2 <- summarise_and_save(
  shrink_contrast(c("Condition", "Resistant_Brequinar", "Resistant_Control")),
  "COMPARISON 2: Brequinar vs Control in Resistant cells",
  "results/DESeq2/DESeq2_Resistant_Brequinar_vs_Control.csv")

comp3 <- summarise_and_save(
  shrink_contrast(c("Condition", "Sensitive_Brequinar", "Sensitive_Control")),
  "COMPARISON 3: Brequinar vs Control in Sensitive cells",
  "results/DESeq2/DESeq2_Sensitive_Brequinar_vs_Control.csv")

# -------------------------------------
# COMPARISON 4: INTERACTION
#   Does brequinar change a gene more, or less, in resistant cells?
# -------------------------------------

dds_interaction <- DESeqDataSetFromMatrix(countData = count_matrix_filtered,
                                          colData   = sample_info,
                                          design    = ~ Resistance + Treatment + Resistance:Treatment)
cat("Running DESeq2 (interaction model)...\n")
dds_interaction <- DESeq(dds_interaction)

cat("\nCoefficients:\n")
print(resultsNames(dds_interaction))
cat("\n")

interaction_term <- "ResistanceResistant.TreatmentBrequinar"
stopifnot(interaction_term %in% resultsNames(dds_interaction))

comp4 <- summarise_and_save(
  lfcShrink(dds_interaction, coef = interaction_term, type = "ashr", quiet = TRUE,
            res = results(dds_interaction, name = interaction_term, alpha = PADJ_CUTOFF)),
  "COMPARISON 4: Interaction (Resistance x Treatment)",
  "results/DESeq2/DESeq2_Interaction_Effect.csv")

cat("Interpretation of the interaction log2FC:\n")
cat("  = (brequinar effect in Resistant) - (brequinar effect in Sensitive)\n")
cat("  Positive: brequinar raises the gene more (or lowers it less) in resistant cells\n")
cat("  Negative: brequinar raises the gene less (or lowers it more) in resistant cells\n\n")

# -------------------------------------
# SUMMARY TABLE
# -------------------------------------

cat("SUMMARY OF ALL COMPARISONS\n")

all_comps <- list(comp1, comp2, comp3, comp4)
summary_stats <- data.frame(
  Comparison = c("Resistant vs Sensitive (untreated)",
                 "Brequinar vs Control (Resistant)",
                 "Brequinar vs Control (Sensitive)",
                 "Interaction (Resistance x Treatment)"),
  Genes_Tested = nrow(count_matrix_filtered),
  do.call(rbind, lapply(all_comps, `[[`, "stats")),
  row.names = NULL
)
summary_stats$Percent_Sig <- round(summary_stats$Sig_padj05 / summary_stats$Genes_Tested * 100, 2)

print(summary_stats)

write.csv(summary_stats, "results/DESeq2/DESeq2_Analysis_Summary.csv", row.names = FALSE)
cat("\nSaved: results/DESeq2/DESeq2_Analysis_Summary.csv\n\n")

# -------------------------------------
# SAVE
#   Same object names that steps 03-07 load.
# -------------------------------------

res_resistant_vs_sensitive_ordered  <- comp1$res
res_resistant_BQ_vs_Control_ordered <- comp2$res
res_sensitive_BQ_vs_Control_ordered <- comp3$res
res_interaction_ordered             <- comp4$res

save(dds, dds_interaction,
     res_resistant_vs_sensitive_ordered,
     res_resistant_BQ_vs_Control_ordered,
     res_sensitive_BQ_vs_Control_ordered,
     res_interaction_ordered,
     sample_info,
     summary_stats,
     file = "results/DESeq2/DESeq2_Results.RData")

cat("Saved: results/DESeq2/DESeq2_Results.RData\n")