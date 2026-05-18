source("scripts/setup.R")

# Loading data
load("results/DESeq2/DESeq2_Results.RData")

# -------------------------------------
# VISUALIZATION OF DEGs
# -------------------------------------

cat("\nCREATING VOLCANO PLOTS\n")

# Function to create volcano plot
create_volcano <- function(res, title, file_name, 
                           pCutoff = 0.05, FCcutoff = 1) {
  
  pdf(file_name, width = 12, height = 10)
  
  p <- EnhancedVolcano(res,
                       lab = rownames(res),
                       x = 'log2FoldChange',
                       y = 'padj',
                       title = title,
                       pCutoff = pCutoff,
                       FCcutoff = FCcutoff,
                       pointSize = 2.0,
                       labSize = 4.0,
                       colAlpha = 0.5,
                       legendPosition = 'right',
                       legendLabSize = 12,
                       legendIconSize = 4.0,
                       drawConnectors = TRUE,
                       widthConnectors = 0.5,
                       max.overlaps = 15)
  
  print(p)
  dev.off()
  
  cat("Saved:", file_name, "\n")
}

cat("Creating volcano plots...\n")

# ---- Volcano plots for all comparisons ----

create_volcano(
  res_resistant_vs_sensitive_ordered,
  "Resistant vs Sensitive (Control)",
  "results/visualization/Volcano_Resistant_vs_Sensitive.pdf"
)

create_volcano(
  res_resistant_BQ_vs_Control_ordered,
  "Brequinar Effect in Resistant Cells",
  "results/visualization/Volcano_Resistant_Brequinar.pdf"
)

create_volcano(
  res_sensitive_BQ_vs_Control_ordered,
  "Brequinar Effect in Sensitive Cells",
  "results/visualization/Volcano_Sensitive_Brequinar.pdf"
)

create_volcano(
  res_interaction_ordered,
  "Interaction: Differential Response to Brequinar",
  "results/visualization/Volcano_Interaction_Effect.pdf"
)

# -------------------------------------
# MA PLOTS 
# -------------------------------------

cat("\nCREATING MA PLOTS\n")

plot_MA_proper <- function(res_obj, title_text) {
  # Extract data
  df <- as.data.frame(res_obj)
  df$significant <- ifelse(!is.na(df$padj) & df$padj < 0.05, "Significant", "Not Significant")
  
  # Remove NA values
  df <- df[!is.na(df$log2FoldChange) & !is.na(df$baseMean), ]
  
  # Create plot
  plot(df$baseMean, 
       df$log2FoldChange,
       log = "x",  # Log scale on x-axis
       pch = 20,
       cex = 0.3,
       col = ifelse(df$significant == "Significant", 
                    rgb(0, 0, 1, 0.5),  # Blue with transparency
                    rgb(0.5, 0.5, 0.5, 0.3)),  # Gray with transparency
       ylim = c(-15, 15),
       xlab = "Mean of normalized counts",
       ylab = "Log2 fold change",
       main = title_text,
       las = 1)
  
  # Add horizontal line at 0
  abline(h = 0, col = "red3", lwd = 2, lty = 2)
  
  # Add lowess smoothing line
  lowess_data <- lowess(log10(df$baseMean + 1), df$log2FoldChange)
  lines(10^lowess_data$x, lowess_data$y, col = "red3", lwd = 2)
  
  # Add legend
  legend("topright", 
         legend = c("Significant (padj < 0.05)", "Not Significant"),
         col = c(rgb(0, 0, 1, 0.5), rgb(0.5, 0.5, 0.5, 0.3)),
         pch = 20,
         pt.cex = 1.5,
         cex = 0.8,
         bg = "white")
  
  # Add counts
  n_sig <- sum(df$significant == "Significant")
  n_total <- nrow(df)
  mtext(paste0(n_sig, " / ", n_total, " genes significant"), 
        side = 3, line = 0.5, cex = 0.8, col = "blue")
}

# Create MA plots
pdf("results/visualization/MA_Plots.pdf", width = 12, height = 10)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

plot_MA_proper(res_resistant_vs_sensitive_ordered,
               "MA Plot: Resistant vs Sensitive")

plot_MA_proper(res_resistant_BQ_vs_Control_ordered,
               "MA Plot: Brequinar in Resistant")

plot_MA_proper(res_sensitive_BQ_vs_Control_ordered,
               "MA Plot: Brequinar in Sensitive")

plot_MA_proper(res_interaction_ordered,
               "MA Plot: Interaction Effect")

par(mfrow = c(1, 1))
dev.off()

cat("Saved: results/visualization/MA_Plots.pdf\n\n")

# -------------------------------------
# HEATMAP OF TOP DEGs
# -------------------------------------

cat("CREATING HEATMAPS OF TOP DEGs\n")

# Function to create heatmap of top DEGs
plot_top_genes_heatmap <- function(res, dds_object, n_genes = 50, title, filename) {
  
  # Get top genes by adjusted p-value
  top_genes <- head(rownames(res[order(res$padj), ]), n_genes)
  
  # Get normalized counts
  norm_counts <- counts(dds_object, normalized = TRUE)
  top_counts <- norm_counts[top_genes, , drop = FALSE]
  
  # Log transform and scale
  log_top_counts <- log2(top_counts + 1)
  scaled_counts <- t(scale(t(log_top_counts)))
  
  # Create annotation
  annotation_col <- data.frame(
    Resistance = sample_info$Resistance,
    Treatment = sample_info$Treatment,
    row.names = colnames(scaled_counts)
  )
  
  annotation_colors <- list(
    Resistance = c(Sensitive = "blue", Resistant = "red"),
    Treatment = c(Control = "gray70", Brequinar = "orange")
  )
  
  # Plot heatmap and save to PDF
  pdf(filename, width = 10, height = 12)
  pheatmap(scaled_counts,
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           show_rownames = TRUE,
           show_colnames = TRUE,
           scale = "none",
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           main = title,
           fontsize = 8,
           fontsize_row = 6,
           color = colorRampPalette(c("blue", "white", "red"))(100))
  dev.off()
  
  cat("Saved:", filename, "\n")
}

cat("Creating heatmaps of top DEGs...\n\n")

# Heatmaps for each comparison
plot_top_genes_heatmap(
  res_resistant_vs_sensitive_ordered,
  dds,  # Using main dds object
  50,
  "Top 50 DEGs: Resistant vs Sensitive",
  "results/visualization/Heatmap_Top50_Resistant_vs_Sensitive.pdf"
)

plot_top_genes_heatmap(
  res_resistant_BQ_vs_Control_ordered,
  dds,  # Using main dds object
  50,
  "Top 50 DEGs: Brequinar in Resistant",
  "results/visualization/Heatmap_Top50_Brequinar_Resistant.pdf"
)

plot_top_genes_heatmap(
  res_sensitive_BQ_vs_Control_ordered,
  dds,  # Use main dds object
  50,
  "Top 50 DEGs: Brequinar in Sensitive",
  "results/visualization/Heatmap_Top50_Brequinar_Sensitive.pdf"
)

plot_top_genes_heatmap(
  res_interaction_ordered,
  dds_interaction,  # Using interaction dds object
  50,
  "Top 50 DEGs: Interaction Effect",
  "results/visualization/Heatmap_Top50_Interaction_Effect.pdf"
)

cat("\nAll heatmaps saved successfully\n\n")