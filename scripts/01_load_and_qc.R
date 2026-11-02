source("scripts/setup.R")

# -------------------------------------
# LOAD COUNT DATA
#   GEO has one count file per sample. I combine them into one matrix after
#   checking that every file lists the same genes in the same order.
# -------------------------------------

cat("\nLOADING AND COMBINING COUNT DATA\n")

count_files <- list.files("data/", pattern = "\\.count\\.out\\.txt$", full.names = TRUE)
stopifnot(length(count_files) == 12)

count_data_list <- list()
gene_names <- NULL

for (file in count_files) {
  sample_name <- sub("\\.count\\.out\\.txt$", "", basename(file))
  sample_name <- sub("^GSM[0-9]+_", "", sample_name)
  counts <- read.table(file, header = TRUE, row.names = 1)

  # Stop if a file's genes are in a different order
  if (is.null(gene_names)) gene_names <- rownames(counts)
  stopifnot(identical(rownames(counts), gene_names))

  count_data_list[[sample_name]] <- counts[, 1]
}

count_matrix <- do.call(cbind, count_data_list)
rownames(count_matrix) <- gene_names

cat("Matrix dimensions:", dim(count_matrix), "\n")
cat("Sample names:", colnames(count_matrix), "\n\n")

# -------------------------------------
# SAMPLE METADATA
#   Groups come from the sample names (WT = sensitive, CR = resistant,
#   C = control, BQ = brequinar), not from file order, so labels can't shift.
#   Sensitive and Control are the reference levels.
# -------------------------------------

sample_info <- data.frame(
  Sample     = colnames(count_matrix),
  Resistance = factor(ifelse(grepl("^OV5CR", colnames(count_matrix)), "Resistant", "Sensitive"),
                      levels = c("Sensitive", "Resistant")),
  Treatment  = factor(ifelse(grepl("BQ", colnames(count_matrix)), "Brequinar", "Control"),
                      levels = c("Control", "Brequinar")),
  row.names  = colnames(count_matrix)
)
sample_info$Condition <- factor(paste(sample_info$Resistance, sample_info$Treatment, sep = "_"))

# 2 x 2 design, 3 replicates per group
stopifnot(all(table(sample_info$Resistance, sample_info$Treatment) == 3))
stopifnot(identical(colnames(count_matrix), rownames(sample_info)))

cat("=== SAMPLE METADATA ===\n")
print(sample_info)
cat("\n")

# -------------------------------------
# QUALITY CONTROL
#   Before testing anything, I want to see that library sizes are similar
#   and that replicates look alike.
# -------------------------------------

cat("QUALITY CONTROL\n")

condition_colors <- c(Sensitive_Control   = "lightblue",
                      Sensitive_Brequinar = "lightgreen",
                      Resistant_Control   = "coral",
                      Resistant_Brequinar = "gold")

library_sizes <- colSums(count_matrix)

pdf("results/QC/QC_Library_Sizes.pdf", width = 10, height = 6)
par(mar = c(8, 4, 4, 2))
barplot(library_sizes / 1e6, las = 2, ylab = "Library size (millions)",
        main = "Library Sizes Across Samples",
        col = condition_colors[as.character(sample_info$Condition)],
        border = NA, cex.names = 0.8)
abline(h = mean(library_sizes / 1e6), col = "red", lty = 2, lwd = 2)
legend("topright",
       legend = c("Sensitive Control", "Sensitive + Brequinar",
                  "Resistant Control", "Resistant + Brequinar", "Mean"),
       fill   = c(unname(condition_colors), NA),
       border = c(rep("black", 4), NA),
       lty    = c(NA, NA, NA, NA, 2), col = c(NA, NA, NA, NA, "red"), cex = 0.7)
invisible(dev.off())

# Keep genes with >= 10 counts in at least 3 samples (the size of one group)
keep <- rowSums(count_matrix >= 10) >= 3
count_matrix_filtered <- count_matrix[keep, ]

cat("=== FILTERING RESULTS ===\n")
cat("Genes before filtering:", nrow(count_matrix), "\n")
cat("Genes after filtering: ", nrow(count_matrix_filtered), "\n")
cat("Genes removed:         ", nrow(count_matrix) - nrow(count_matrix_filtered), "\n\n")

# log-CPM is only for the QC plots; DESeq2 uses the raw counts
dge <- DGEList(counts = count_matrix_filtered)
dge <- calcNormFactors(dge, method = "TMM")
cpm_data <- cpm(dge, log = TRUE, prior.count = 1)

cat("=== DATA QUALITY CHECKS ===\n")
cat("Any NA in log-CPM? ", any(is.na(cpm_data)), "\n")
cat("Any Inf in log-CPM?", any(is.infinite(cpm_data)), "\n")
cat("log-CPM range:     ", round(range(cpm_data), 2), "\n\n")

# -------------------------------------
# PCA
#   Do samples split by resistance and treatment, and do replicates stay together?
# -------------------------------------

pca <- prcomp(t(cpm_data))
pca_df <- cbind(as.data.frame(pca$x), sample_info)
variance_explained <- summary(pca)$importance[2, ]

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Resistance, shape = Treatment)) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = Sample), size = 3, max.overlaps = Inf,
                  box.padding = 0.5, show.legend = FALSE) +
  theme_bw() +
  labs(title = "PCA of Samples (TMM log-CPM)",
       x = paste0("PC1 (", round(variance_explained[1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2] * 100, 1), "%)")) +
  scale_color_manual(values = c(Sensitive = "blue", Resistant = "red"))

pdf("results/QC/QC_PCA_Plot.pdf", width = 10, height = 8)
print(p)
invisible(dev.off())

# PNG copy for the README
ggsave("results/figures/pca_plot.png", p, width = 10, height = 8, dpi = 150)

# -------------------------------------
# SAMPLE DISTANCE & CORRELATION HEATMAPS
#   Same question as the PCA, from a different angle.
# -------------------------------------

sample_dist <- dist(t(cpm_data))
annotation_col <- data.frame(Resistance = sample_info$Resistance,
                             Treatment  = sample_info$Treatment,
                             row.names  = rownames(sample_info))
annotation_colors <- list(Resistance = c(Sensitive = "blue", Resistant = "red"),
                          Treatment  = c(Control = "gray70", Brequinar = "orange"))

pdf("results/QC/QC_Sample_Distance_Heatmap.pdf", width = 10, height = 9)
pheatmap(as.matrix(sample_dist),
         clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         main = "Sample-to-Sample Distance (TMM log-CPM)",
         fontsize = 10)
invisible(dev.off())

pdf("results/QC/QC_Sample_Correlation_Heatmap.pdf", width = 10, height = 9)
pheatmap(cor(cpm_data),
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         main = "Sample Correlation (TMM log-CPM)",
         display_numbers = TRUE, number_format = "%.2f",
         fontsize = 10, fontsize_number = 8)
invisible(dev.off())

# -------------------------------------
# SAVE
# -------------------------------------

save(count_matrix_filtered, sample_info, cpm_data, dge,
     file = "results/QC/Filtered_Data.RData")

cat("Saved: results/QC/Filtered_Data.RData\n")