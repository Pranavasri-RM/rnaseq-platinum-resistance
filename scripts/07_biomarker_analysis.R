source("scripts/setup.R")

load("results/DESeq2/DESeq2_Results.RData")
load("results/DESeq2/DESeq2_Annotated_Results.RData")
load("results/enrichment/Pyrimidine_Metabolism_Results.RData")

# ============================================
# CANDIDATE RESISTANCE MARKERS
#   Genes that strongly separate untreated resistant and sensitive cells:
#   padj < 0.01, |log2FC| > 2, baseMean > 100. I call them candidates because
#   they come from one cell-line pair and would need checking in other models.
# ============================================

PADJ_MARKER <- 0.01
LFC_MARKER  <- 2
MIN_MEAN    <- 100
OUT_DIR     <- "results/biomarkers"

cat("\nCANDIDATE RESISTANCE MARKERS (untreated Resistant vs Sensitive)\n")

markers <- res_resistant_annotated %>%
  filter(!is.na(padj), padj < PADJ_MARKER,
         abs(log2FoldChange) > LFC_MARKER,
         baseMean > MIN_MEAN) %>%
  arrange(padj)

# Check 1: do the 3 resistant and 3 sensitive replicates separate completely?
norm_counts <- counts(dds, normalized = TRUE)
res_ctrl <- rownames(sample_info)[sample_info$Condition == "Resistant_Control"]
sen_ctrl <- rownames(sample_info)[sample_info$Condition == "Sensitive_Control"]

markers$groups_separate <- sapply(seq_len(nrow(markers)), function(i) {
  g <- markers$gene_symbol[i]
  r <- norm_counts[g, res_ctrl]; s <- norm_counts[g, sen_ctrl]
  if (markers$log2FoldChange[i] > 0) min(r) > max(s) else max(r) < min(s)
})

# Check 2: is the difference still there when both lines are treated with brequinar?
res_bq_pair <- lfcShrink(dds,
                         contrast = c("Condition", "Resistant_Brequinar", "Sensitive_Brequinar"),
                         res = results(dds, contrast = c("Condition", "Resistant_Brequinar",
                                                         "Sensitive_Brequinar")),
                         type = "ashr", quiet = TRUE)
markers$log2FC_under_BQ <- res_bq_pair[markers$gene_symbol, "log2FoldChange"]
markers$padj_under_BQ   <- res_bq_pair[markers$gene_symbol, "padj"]
markers$persists_under_BQ <- !is.na(markers$padj_under_BQ) &
  markers$padj_under_BQ < 0.05 &
  sign(markers$log2FC_under_BQ) == sign(markers$log2FoldChange)

markers <- markers %>%
  relocate(gene_symbol, log2FoldChange, padj, baseMean,
           groups_separate, log2FC_under_BQ, padj_under_BQ, persists_under_BQ)

cat("Candidates:", nrow(markers),
    "(", sum(markers$log2FoldChange > 0), "up,", sum(markers$log2FoldChange < 0), "down in resistant )\n")
cat("  No overlap between resistant and sensitive replicates:", sum(markers$groups_separate), "\n")
cat("  Same direction and padj < 0.05 under brequinar too:   ", sum(markers$persists_under_BQ), "\n\n")

cat("Top 15 by padj:\n")
print(head(markers %>%
             transmute(gene_symbol,
                       log2FC = round(log2FoldChange, 2),
                       padj = signif(padj, 2),
                       baseMean = round(baseMean),
                       groups_separate, persists_under_BQ), 15),
      row.names = FALSE)
cat("\n")

write.csv(markers, file.path(OUT_DIR, "Resistance_Marker_Candidates.csv"), row.names = FALSE)
cat("Saved:", file.path(OUT_DIR, "Resistance_Marker_Candidates.csv"), "\n")

# ---- Expression of the top 12 up and top 12 down candidates in all four groups ----
vsd <- vst(dds, blind = FALSE)

top_up   <- head(markers$gene_symbol[markers$log2FoldChange > 0], 12)
top_down <- head(markers$gene_symbol[markers$log2FoldChange < 0], 12)

plot_markers <- function(genes, title) {
  if (length(genes) == 0) return(NULL)
  df <- as.data.frame(t(assay(vsd)[genes, , drop = FALSE]))
  df$Sample <- rownames(df)
  df <- tidyr::pivot_longer(df, -Sample, names_to = "Gene", values_to = "VST")
  df <- left_join(df, sample_info %>% select(Sample, Resistance, Treatment), by = "Sample")
  df$Gene <- factor(df$Gene, levels = genes)

  ggplot(df, aes(x = Resistance, y = VST, colour = Treatment)) +
    geom_point(position = position_dodge(width = 0.5), size = 2) +
    stat_summary(aes(group = Treatment), fun = mean, geom = "crossbar",
                 width = 0.35, position = position_dodge(width = 0.5)) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
    scale_colour_manual(values = c(Control = "grey40", Brequinar = "darkorange")) +
    theme_bw(base_size = 10) +
    labs(title = title, y = "Expression (VST)", x = NULL)
}

pdf(file.path(OUT_DIR, "Resistance_Marker_Candidates_Expression.pdf"), width = 11, height = 8.5)
p_up   <- plot_markers(top_up,   "Top candidates higher in resistant cells")
p_down <- plot_markers(top_down, "Top candidates lower in resistant cells")
if (!is.null(p_up))   print(p_up)
if (!is.null(p_down)) print(p_down)
invisible(dev.off())
cat("Saved:", file.path(OUT_DIR, "Resistance_Marker_Candidates_Expression.pdf"), "\n")

save(markers, res_bq_pair, file = file.path(OUT_DIR, "Resistance_Marker_Results.RData"))

# ============================================
# PIPELINE SUMMARY
#   Headline numbers from each step in one table.
# ============================================

cat("\nFINAL PIPELINE SUMMARY\n")

final_summary <- data.frame(
  Metric = c(
    "Genes tested (after low-count filter)",
    "DE genes, Resistant vs Sensitive (padj < 0.05)",
    "DE genes, Brequinar in Resistant (padj < 0.05)",
    "DE genes, Brequinar in Sensitive (padj < 0.05)",
    "Interaction genes (padj < 0.05)",
    "Candidate resistance markers (padj < 0.01, |log2FC| > 2, baseMean > 100)",
    "  of which replicates do not overlap",
    "  of which persist under brequinar",
    "Curated pyrimidine genes tested"
  ),
  Value = c(
    nrow(res_resistant_annotated),
    summary_stats$Sig_padj05,
    nrow(markers),
    sum(markers$groups_separate),
    sum(markers$persists_under_BQ),
    nrow(pyrimidine_wide)
  )
)
print(final_summary, row.names = FALSE)

write.csv(final_summary, file.path(OUT_DIR, "Final_Pipeline_Summary.csv"), row.names = FALSE)
cat("\nSaved:", file.path(OUT_DIR, "Final_Pipeline_Summary.csv"), "\n")

# Save package versions so the run can be reproduced
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
cat("Saved: results/sessionInfo.txt\n")