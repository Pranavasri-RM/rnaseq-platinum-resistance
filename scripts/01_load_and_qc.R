source("scripts/setup.R")

# Input files
count_files <- list.files("data/", pattern = "count.out.txt$", full.names = TRUE)

# -------------------------------------
# LOAD AND COMBINE COUNT DATA
# -------------------------------------

cat("\nLOADING AND COMBINING COUNT DATA\n")

count_data_list <- list()
first_file <- read.table(count_files[1], header = TRUE, row.names = 1)
gene_names <- rownames(first_file)

for (file in count_files) {
  sample_name <- gsub("GSM[0-9]+_", "", file)
  sample_name <- gsub(".count.out.txt", "", sample_name)
  counts <- read.table(file, header = TRUE, row.names = 1)
  count_data_list[[sample_name]] <- counts[, 1]
}

count_matrix <- do.call(cbind, count_data_list)
rownames(count_matrix) <- gene_names
colnames(count_matrix) <- names(count_data_list)

cat("First 20 gene names:\n")
print(rownames(count_matrix)[1:20])
cat("\nMatrix dimensions:", dim(count_matrix), "\n")
cat("Sample names:", colnames(count_matrix), "\n")
head(count_matrix)

# -------------------------------------
#   CREATE SAMPLE METADATA
# -------------------------------------

sample_info <- data.frame(
  Sample = colnames(count_matrix),
  Resistance = c(rep("Sensitive", 6), rep("Resistant", 6)),
  Treatment = rep(c(rep("Control", 3), rep("Brequinar", 3)), 2),
  row.names = colnames(count_matrix)
)

sample_info$Condition <- factor(paste(sample_info$Resistance, 
                                      sample_info$Treatment, 
                                      sep = "_"))

cat("=== SAMPLE METADATA ===\n")
print(sample_info)
cat("\nCount matrix columns match metadata rows:", 
    all(colnames(count_matrix) == rownames(sample_info)), "\n\n")

# -------------------------------------
#   QUALITY CONTROL
# -------------------------------------

cat("QUALITY CONTROL\n")

library_sizes <- colSums(count_matrix)

pdf("results/QC/QC_Library_Sizes.pdf", width = 10, height = 6)
barplot(library_sizes/1e6, las = 2, ylab = "Library Size (millions)",
        main = "Library Sizes Across Samples",
        col = c(rep("lightblue", 3), rep("lightgreen", 3), rep("coral", 3), rep("gold", 3)),
        border = NA, cex.names = 0.8)
abline(h = mean(library_sizes/1e6), col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("Sensitive Control", "Sensitive + Brequinar", 
                              "Resistant Control", "Resistant + Brequinar", "Mean"),
       fill = c("lightblue", "lightgreen", "coral", "gold", NA),
       border = c("black", "black", "black", "black", NA),
       lty = c(NA, NA, NA, NA, 2), col = c(NA, NA, NA, NA, "red"), cex = 0.7)
dev.off()

# Filtering low-count genes
keep <- rowSums(count_matrix >= 10) >= 3
count_matrix_filtered <- count_matrix[keep, ]

cat("=== FILTERING RESULTS ===\n")
cat("Genes before filtering:", nrow(count_matrix), "\n")
cat("Genes after filtering:", nrow(count_matrix_filtered), "\n")
cat("Genes removed:", nrow(count_matrix) - nrow(count_matrix_filtered), "\n\n")

# Normalization (CPM) for QC
dge <- DGEList(counts = count_matrix_filtered)
dge <- calcNormFactors(dge, method = "TMM")
cpm_data <- cpm(dge, log = TRUE, prior.count = 1)

cat("=== DATA QUALITY CHECKS ===\n")
cat("Any NA in CPM?", any(is.na(cpm_data)), "\n")
cat("Any Inf in CPM?", any(is.infinite(cpm_data)), "\n")
cat("CPM range:", round(range(cpm_data), 2), "\n\n")

# PCA
pca <- prcomp(t(cpm_data))
pca_df <- as.data.frame(pca$x)
pca_df <- cbind(pca_df, sample_info)
variance_explained <- summary(pca)$importance[2,]

pdf("results/QC/QC_PCA_Plot.pdf", width = 10, height = 8)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Resistance, shape = Treatment)) +
  geom_point(size = 5) +
  geom_text_repel(aes(label = Sample), size = 3, 
                  max.overlaps = Inf,  # Show all labels
                  box.padding = 0.5) +
  theme_bw() +
  labs(title = "PCA of Sample Groups",
       x = paste0("PC1 (", round(variance_explained[1]*100, 1), "%)"),
       y = paste0("PC2 (", round(variance_explained[2]*100, 1), "%)")) +
  scale_color_manual(values = c("Sensitive" = "blue", "Resistant" = "red"))
dev.off()

# Distance & Correlation Heatmaps
sample_dist <- dist(t(cpm_data))
sample_dist_matrix <- as.matrix(sample_dist)
annotation_col <- data.frame(Resistance = sample_info$Resistance,
                             Treatment = sample_info$Treatment,
                             row.names = rownames(sample_info))
annotation_colors <- list(Resistance = c(Sensitive="blue", Resistant="red"),
                          Treatment = c(Control="gray70", Brequinar="orange"))

pdf("results/QC/QC_Sample_Distance_Heatmap.pdf", width = 10, height = 9)
pheatmap(sample_dist_matrix, clustering_distance_rows = sample_dist,
         clustering_distance_cols = sample_dist,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         main = "Sample-to-Sample Distance Heatmap",
         fontsize = 10)
dev.off()

cor_matrix <- cor(cpm_data)
pdf("results/QC/QC_Sample_Correlation_Heatmap.pdf", width = 10, height = 9)
pheatmap(cor_matrix, annotation_col = annotation_col, 
         annotation_colors = annotation_colors,
         main = "Sample Correlation Heatmap",
         display_numbers = TRUE, number_format = "%.2f", fontsize = 10, fontsize_number = 8)
dev.off()

# Save filtered objects
save(count_matrix_filtered, sample_info, cpm_data, dge,
     file = "results/QC/Filtered_Data.RData")