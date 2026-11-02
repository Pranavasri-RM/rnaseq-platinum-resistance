source("scripts/setup.R")

load("results/DESeq2/DESeq2_Results.RData")

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1

# -------------------------------------
# VOLCANO PLOTS
#   y-axis is padj, since that's what the DE calls use.
#   Genes with padj = 0 get pinned to the top (EnhancedVolcano warns about this).
# -------------------------------------

cat("\nCREATING VOLCANO PLOTS\n")

make_volcano <- function(res, title) {
  EnhancedVolcano(res,
                  lab = rownames(res),
                  x = "log2FoldChange",
                  y = "padj",
                  title = title,
                  subtitle = paste0("padj < ", PADJ_CUTOFF, ", |log2FC| > ", LFC_CUTOFF),
                  ylab = bquote(~-Log[10] ~ "adjusted" ~ italic(P)),
                  legendLabels = c("NS", "Log2 FC", "padj", "padj & Log2 FC"),
                  pCutoff = PADJ_CUTOFF,
                  FCcutoff = LFC_CUTOFF,
                  pointSize = 2.0,
                  labSize = 4.0,
                  colAlpha = 0.5,
                  legendPosition = "right",
                  legendLabSize = 12,
                  legendIconSize = 4.0,
                  drawConnectors = TRUE,
                  widthConnectors = 0.5,
                  max.overlaps = 15)
}

save_volcano <- function(res, title, pdf_file, png_file = NULL) {
  p <- make_volcano(res, title)
  pdf(pdf_file, width = 12, height = 10)
  print(p)
  invisible(dev.off())
  cat("Saved:", pdf_file, "\n")
  if (!is.null(png_file)) {
    ggsave(png_file, p, width = 12, height = 10, dpi = 150)
    cat("Saved:", png_file, "\n")
  }
}

save_volcano(res_resistant_vs_sensitive_ordered,
             "Resistant vs Sensitive (untreated)",
             "results/visualization/Volcano_Resistant_vs_Sensitive.pdf",
             "results/figures/volcano_resistant_vs_sensitive.png")   # used in README

save_volcano(res_resistant_BQ_vs_Control_ordered,
             "Brequinar vs Control in Resistant Cells",
             "results/visualization/Volcano_Resistant_Brequinar.pdf")

save_volcano(res_sensitive_BQ_vs_Control_ordered,
             "Brequinar vs Control in Sensitive Cells",
             "results/visualization/Volcano_Sensitive_Brequinar.pdf")

save_volcano(res_interaction_ordered,
             "Interaction: Brequinar Effect in Resistant minus Sensitive",
             "results/visualization/Volcano_Interaction_Effect.pdf")

# -------------------------------------
# MA PLOTS
#   Mainly to check for fold-change bias at low counts.
#   The red lowess line should stay close to zero.
# -------------------------------------

cat("\nCREATING MA PLOTS\n")

plot_ma <- function(res, title_text, ylim = c(-10, 10)) {
  df <- as.data.frame(res)
  df <- df[!is.na(df$log2FoldChange) & df$baseMean > 0, ]
  sig <- !is.na(df$padj) & df$padj < PADJ_CUTOFF

  # Points past the y-limits are drawn at the edge as triangles
  y <- pmax(pmin(df$log2FoldChange, ylim[2]), ylim[1])
  clipped <- y != df$log2FoldChange

  plot(df$baseMean, y, log = "x",
       pch = ifelse(clipped, 17, 20), cex = ifelse(clipped, 0.6, 0.3),
       col = ifelse(sig, rgb(0, 0, 1, 0.5), rgb(0.5, 0.5, 0.5, 0.3)),
       ylim = ylim, las = 1,
       xlab = "Mean of normalized counts", ylab = "Log2 fold change",
       main = title_text)
  abline(h = 0, col = "red3", lwd = 2, lty = 2)

  lw <- lowess(log10(df$baseMean), df$log2FoldChange)
  lines(10^lw$x, lw$y, col = "red3", lwd = 2)

  legend("topright",
         legend = c(paste0("padj < ", PADJ_CUTOFF), "Not significant", "Beyond y-limit"),
         col = c(rgb(0, 0, 1, 0.5), rgb(0.5, 0.5, 0.5, 0.3), "black"),
         pch = c(20, 20, 17), pt.cex = 1.2, cex = 0.8, bg = "white")
  mtext(paste0(sum(sig), " / ", nrow(df), " genes significant"),
        side = 3, line = 0.5, cex = 0.8, col = "blue")
}

pdf("results/visualization/MA_Plots.pdf", width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(4, 4, 4, 1))
plot_ma(res_resistant_vs_sensitive_ordered,  "Resistant vs Sensitive (untreated)")
plot_ma(res_resistant_BQ_vs_Control_ordered, "Brequinar in Resistant")
plot_ma(res_sensitive_BQ_vs_Control_ordered, "Brequinar in Sensitive")
plot_ma(res_interaction_ordered,             "Interaction (Resistance x Treatment)")
invisible(dev.off())
cat("Saved: results/visualization/MA_Plots.pdf\n\n")

# -------------------------------------
# HEATMAPS OF TOP 50 GENES
#   VST counts, z-scored per gene and capped at +/-3
#   so a few extreme genes don't flatten the colour scale.
# -------------------------------------

cat("CREATING HEATMAPS OF TOP DEGs\n")

vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

annotation_col <- data.frame(Resistance = sample_info$Resistance,
                             Treatment  = sample_info$Treatment,
                             row.names  = rownames(sample_info))
annotation_colors <- list(Resistance = c(Sensitive = "blue", Resistant = "red"),
                          Treatment  = c(Control = "gray70", Brequinar = "orange"))
stopifnot(identical(colnames(vst_mat), rownames(annotation_col)))

plot_top_heatmap <- function(res, title, filename, n_genes = 50) {
  res_sig <- res[!is.na(res$padj) & res$padj < PADJ_CUTOFF, ]
  top_genes <- head(rownames(res_sig[order(res_sig$padj), ]), n_genes)
  if (length(top_genes) < 2) {
    cat("Fewer than 2 significant genes, skipping:", filename, "\n")
    return(invisible(NULL))
  }

  z <- t(scale(t(vst_mat[top_genes, , drop = FALSE])))
  z[z >  3] <-  3     # cap extreme z-scores so the color scale stays readable
  z[z < -3] <- -3

  pdf(filename, width = 10, height = 12)
  pheatmap(z,
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           cluster_rows = TRUE, cluster_cols = TRUE,
           main = paste0(title, " (VST, row z-score)"),
           fontsize = 8, fontsize_row = 6,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(-3, 3, length.out = 101))
  invisible(dev.off())
  cat("Saved:", filename, "\n")
}

plot_top_heatmap(res_resistant_vs_sensitive_ordered,  "Top 50: Resistant vs Sensitive",
                 "results/visualization/Heatmap_Top50_Resistant_vs_Sensitive.pdf")
plot_top_heatmap(res_resistant_BQ_vs_Control_ordered, "Top 50: Brequinar in Resistant",
                 "results/visualization/Heatmap_Top50_Brequinar_Resistant.pdf")
plot_top_heatmap(res_sensitive_BQ_vs_Control_ordered, "Top 50: Brequinar in Sensitive",
                 "results/visualization/Heatmap_Top50_Brequinar_Sensitive.pdf")
plot_top_heatmap(res_interaction_ordered,             "Top 50: Interaction",
                 "results/visualization/Heatmap_Top50_Interaction_Effect.pdf")

cat("\nAll visualizations saved\n")