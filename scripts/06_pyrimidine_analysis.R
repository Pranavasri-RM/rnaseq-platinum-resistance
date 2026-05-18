source("scripts/setup.R")

# From DESeq step
load("results/DESeq2/DESeq2_Results.RData")

# From Annotation step
load("results/DESeq2/DESeq2_Annotated_Results.RData")

# From pyrimidine step
load("results/enrichment/Pathway_Enrichment_Results.RData")

# ============================================
#  PYRIMIDINE METABOLISM
# ============================================

cat("\nTARGETED ANALYSIS: PYRIMIDINE METABOLISM\n")
cat("Key pathway for understanding Brequinar's mechanism\n")

# Define pyrimidine metabolism genes (key enzymes and regulators)
pyrimidine_keywords <- c("DHODH",   # Dihydroorotate dehydrogenase (direct target)
                         "CAD",     # Carbamoyl-phosphate synthetase
                         "UMPS",    # Uridine 5'-monophosphate synthase
                         "CTPS",    # CTP synthase
                         "TYMS",    # Thymidylate synthase
                         "UCK",     # Uridine-cytidine kinase
                         "CMPK",    # Cytidine monophosphate kinase
                         "NME",     # Nucleoside diphosphate kinase
                         "RRM",     # Ribonucleotide reductase
                         "TK",      # Thymidine kinase
                         "POLD",    # DNA polymerase delta
                         "PNP",     # Purine nucleoside phosphorylase
                         "NT5E")    # 5'-nucleotidase

cat("Pyrimidine pathway genes of interest:\n")
cat(paste(pyrimidine_keywords, collapse = ", "), "\n\n")

# Function to extract and analyze pyrimidine metabolism genes
extract_pyrimidine_genes <- function(res_annotated, comparison_name) {
  
  cat("Searching for pyrimidine metabolism genes in:", comparison_name, "\n")
  
  # Extracting genes matching pyrimidine keywords
  pyrimidine_genes <- res_annotated[
    grepl(paste(pyrimidine_keywords, collapse = "|"), 
          res_annotated$gene_symbol,
          ignore.case = TRUE), ]
  
  n_found <- nrow(pyrimidine_genes)
  cat("  Found:", n_found, "genes\n")
  
  if (n_found == 0) {
    cat("No pyrimidine genes found\n\n")
    return(NULL)
  }
  
  # Sorting by adjusted p-value
  pyrimidine_genes <- pyrimidine_genes[order(pyrimidine_genes$padj), ]
  
  # Counting significant genes
  n_sig <- sum(pyrimidine_genes$padj < 0.05, na.rm = TRUE)
  cat("  Significant (padj < 0.05):", n_sig, "\n\n")
  
  return(pyrimidine_genes)
}

# ============================================
# EXTRACT PYRIMIDINE GENES FOR ALL COMPARISONS
# ============================================

cat("EXTRACTING PYRIMIDINE GENES:\n")

pyrimidine_resistant <- extract_pyrimidine_genes(
  res_resistant_annotated,
  "Resistant vs Sensitive (Control)"
)

pyrimidine_BQ_resistant <- extract_pyrimidine_genes(
  res_BQ_resistant_annotated,
  "Brequinar in Resistant Cells"
)

pyrimidine_BQ_sensitive <- extract_pyrimidine_genes(
  res_BQ_sensitive_annotated,
  "Brequinar in Sensitive Cells"
)

pyrimidine_interaction <- extract_pyrimidine_genes(
  res_interaction_annotated,
  "Interaction Effect"
)

# ============================================
# DISPLAY RESULTS
# ============================================

cat("PYRIMIDINE METABOLISM GENES - RESULTS\n")

display_pyrimidine_results <- function(pyrimidine_df, comparison_name) {
  if (is.null(pyrimidine_df)) {
    cat("---", comparison_name, "---\n")
    cat("No pyrimidine genes found\n\n")
    return()
  }
  
  cat("---", comparison_name, "---\n")
  cat("Top genes by adjusted p-value:\n\n")
  
  # Display key columns
  display_df <- pyrimidine_df[, c("gene_symbol", "baseMean", "log2FoldChange", 
                                  "padj", "gene_biotype")]
  
  # Rename for clarity
  colnames(display_df) <- c("Gene", "Base Mean", "log2FC", "Padj", "Biotype")
  
  print(display_df)
  cat("\n")
}

display_pyrimidine_results(pyrimidine_resistant, 
                           "Resistant vs Sensitive (Control)")
display_pyrimidine_results(pyrimidine_BQ_resistant,
                           "Brequinar in Resistant Cells")
display_pyrimidine_results(pyrimidine_BQ_sensitive,
                           "Brequinar in Sensitive Cells")
display_pyrimidine_results(pyrimidine_interaction,
                           "Interaction Effect")

# ============================================
# SAVE PYRIMIDINE RESULTS
# ============================================

cat("SAVING PYRIMIDINE METABOLISM RESULTS\n")

if (!is.null(pyrimidine_resistant)) {
  write.csv(pyrimidine_resistant, 
            "results/enrichment/Pyrimidine_Metabolism_Resistant_vs_Sensitive.csv",
            row.names = FALSE)
  cat("Saved: Pyrimidine_Metabolism_Resistant_vs_Sensitive.csv\n")
}

if (!is.null(pyrimidine_BQ_resistant)) {
  write.csv(pyrimidine_BQ_resistant,
            "results/enrichment/Pyrimidine_Metabolism_Brequinar_Resistant.csv",
            row.names = FALSE)
  cat("Saved: Pyrimidine_Metabolism_Brequinar_Resistant.csv\n")
}

if (!is.null(pyrimidine_BQ_sensitive)) {
  write.csv(pyrimidine_BQ_sensitive,
            "results/enrichment/Pyrimidine_Metabolism_Brequinar_Sensitive.csv",
            row.names = FALSE)
  cat("Saved: Pyrimidine_Metabolism_Brequinar_Sensitive.csv\n")
}

if (!is.null(pyrimidine_interaction)) {
  write.csv(pyrimidine_interaction,
            "results/enrichment/Pyrimidine_Metabolism_Interaction_Effect.csv",
            row.names = FALSE)
  cat("Saved: Pyrimidine_Metabolism_Interaction_Effect.csv\n")
}

# ============================================
# VISUALIZATION: HEATMAP OF PYRIMIDINE GENES
# ============================================

cat("CREATING HEATMAP VISUALIZATIONS\n")

# Function to create heatmap for pyrimidine genes
plot_pyrimidine_heatmap <- function(res_df, dds_obj, comparison_name, filename) {
  
  if (is.null(res_df) || nrow(res_df) == 0) {
    cat("No genes to visualize for:", comparison_name, "\n")
    return()
  }
  
  pyrimidine_genes <- res_df$gene_symbol
  
  tryCatch({
    # Get normalized counts
    norm_counts <- counts(dds_obj, normalized = TRUE)
    pyr_counts <- norm_counts[pyrimidine_genes, , drop = FALSE]
    
    # Check if any genes found in counts
    if (nrow(pyr_counts) == 0) {
      cat("Genes not found in count matrix for:", comparison_name, "\n")
      return()
    }
    
    # Log transform and scale
    log_pyr_counts <- log2(pyr_counts + 1)
    scaled_counts <- t(scale(t(log_pyr_counts)))
    
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
    
    # Create heatmap with adjusted dimensions
    pdf(filename, width = 12, height = 10)  # Increased width and height
    
    pheatmap(scaled_counts,
             annotation_col = annotation_col,
             annotation_colors = annotation_colors,
             show_rownames = TRUE,
             show_colnames = TRUE,
             scale = "none",
             cluster_cols = TRUE,
             cluster_rows = TRUE,
             main = paste("Pyrimidine Metabolism Genes:\n", comparison_name),
             fontsize = 10,
             fontsize_row = 5,      # Smaller row names
             fontsize_col = 6,      # Column name size
             cellwidth = 25,        # Wider cells
             cellheight = 4,        # Taller cells
             treeheight_row = 40,   # Height of row dendrogram
             treeheight_col = 40,   # Height of column dendrogram
             angle_col = 45,        # Angle column names at 45 degrees
             color = colorRampPalette(c("blue", "white", "red"))(100))
    
    dev.off()
    cat("Saved:", filename, "\n")
    
  }, error = function(e) {
    cat("Error creating heatmap:", conditionMessage(e), "\n")
  })
}

cat("Creating heatmaps for pyrimidine genes...\n\n")

plot_pyrimidine_heatmap(pyrimidine_resistant, dds,
                        "Resistant vs Sensitive (Control)",
                        "results/enrichment/Heatmap_Pyrimidine_Resistant_vs_Sensitive.pdf")

plot_pyrimidine_heatmap(pyrimidine_BQ_resistant, dds,
                        "Brequinar in Resistant Cells",
                        "results/enrichment/Heatmap_Pyrimidine_Brequinar_Resistant.pdf")

plot_pyrimidine_heatmap(pyrimidine_BQ_sensitive, dds,
                        "Brequinar in Sensitive Cells",
                        "results/enrichment/Heatmap_Pyrimidine_Brequinar_Sensitive.pdf")

plot_pyrimidine_heatmap(pyrimidine_interaction, dds_interaction,
                        "Interaction Effect",
                        "results/enrichment/Heatmap_Pyrimidine_Interaction_Effect.pdf")

# ============================================
# SUMMARY STATISTICS
# ============================================

cat("PYRIMIDINE METABOLISM - SUMMARY\n")

summary_pyrimidine <- data.frame(
  Comparison = c("Resistant vs Sensitive",
                 "Brequinar in Resistant",
                 "Brequinar in Sensitive",
                 "Interaction Effect"),
  Genes_Found = c(
    ifelse(is.null(pyrimidine_resistant), 0, nrow(pyrimidine_resistant)),
    ifelse(is.null(pyrimidine_BQ_resistant), 0, nrow(pyrimidine_BQ_resistant)),
    ifelse(is.null(pyrimidine_BQ_sensitive), 0, nrow(pyrimidine_BQ_sensitive)),
    ifelse(is.null(pyrimidine_interaction), 0, nrow(pyrimidine_interaction))
  ),
  Significant = c(
    ifelse(is.null(pyrimidine_resistant), 0, 
           sum(pyrimidine_resistant$padj < 0.05, na.rm = TRUE)),
    ifelse(is.null(pyrimidine_BQ_resistant), 0,
           sum(pyrimidine_BQ_resistant$padj < 0.05, na.rm = TRUE)),
    ifelse(is.null(pyrimidine_BQ_sensitive), 0,
           sum(pyrimidine_BQ_sensitive$padj < 0.05, na.rm = TRUE)),
    ifelse(is.null(pyrimidine_interaction), 0,
           sum(pyrimidine_interaction$padj < 0.05, na.rm = TRUE))
  )
)

print(summary_pyrimidine)

write.csv(summary_pyrimidine,
          "results/enrichment/Pyrimidine_Metabolism_Summary.csv",
          row.names = FALSE)

cat("\nSaved: results/enrichment/Pyrimidine_Metabolism_Summary.csv\n\n")

# ============================================
# SAVE PYRIMIDINE OBJECTS
# ============================================

save(pyrimidine_resistant,
     pyrimidine_BQ_resistant,
     pyrimidine_BQ_sensitive,
     pyrimidine_interaction,
     summary_pyrimidine,
     file = "results/enrichment/Pyrimidine_Metabolism_Results.RData")

cat("All pyrimidine results saved: results/enrichment/Pyrimidine_Metabolism_Results.RData\n")
