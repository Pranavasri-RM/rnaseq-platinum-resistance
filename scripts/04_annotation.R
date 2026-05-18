source("scripts/setup.R")

# Loading data
load("results/DESeq2/DESeq2_Results.RData")

# ============================================
# GENE ANNOTATION AND ID CONVERSION
# ============================================

cat("\nGENE ANNOTATION - GENE SYMBOLS TO ENSEMBL IDs\n")

# Function to add gene annotations from biomaRt
add_gene_annotations <- function(res_df, comparison_name) {
  cat("Annotating genes for:", comparison_name, "\n")
  
  # Get gene symbols from rownames (your data appears to already have gene symbols)
  gene_symbols <- rownames(res_df)
  
  tryCatch({
    # Connect to biomaRt
    cat("  Connecting to Ensembl database...\n")
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    # Retrieve gene annotations using gene SYMBOLS as the filter
    cat("Retrieving annotations...\n")
    annotations <- getBM(
      attributes = c("ensembl_gene_id", 
                     "external_gene_name", 
                     "description",
                     "gene_biotype",
                     "chromosome_name",
                     "start_position",
                     "end_position"),
      filters = "external_gene_name",  # Changed from ensembl_gene_id
      values = gene_symbols,           # Using gene symbols
      mart = mart
    )
    
    # Convert results to data frame
    res_df_converted <- as.data.frame(res_df)
    res_df_converted$gene_symbol <- gene_symbols
    
    # Merge annotations with results
    res_annotated <- merge(res_df_converted, 
                           annotations, 
                           by.x = "gene_symbol",
                           by.y = "external_gene_name",
                           all.x = TRUE)
    
    # Sort by adjusted p-value (maintain order)
    res_annotated <- res_annotated[order(res_annotated$padj, na.last = NA), ]
    
    cat("  Successfully annotated genes\n")
    cat("  Genes with Ensembl IDs:", sum(!is.na(res_annotated$ensembl_gene_id)), "\n")
    cat("  Genes matched:", nrow(annotations), "\n\n")
    
    return(res_annotated)
    
  }, error = function(e) {
    cat("Error in annotation:", conditionMessage(e), "\n")
    cat("Returning results without annotations\n\n")
    res_df_converted <- as.data.frame(res_df)
    res_df_converted$gene_symbol <- gene_symbols
    return(res_df_converted)
  })
}

# ---- Annotate all comparison results ----

cat("Processing Comparison 1: Resistant vs Sensitive\n")
res_resistant_annotated <- add_gene_annotations(
  res_resistant_vs_sensitive_ordered,
  "Resistant vs Sensitive (Control)"
)

cat("Processing Comparison 2: Brequinar in Resistant Cells\n")
res_BQ_resistant_annotated <- add_gene_annotations(
  res_resistant_BQ_vs_Control_ordered,
  "Brequinar in Resistant Cells"
)

cat("Processing Comparison 3: Brequinar in Sensitive Cells\n")
res_BQ_sensitive_annotated <- add_gene_annotations(
  res_sensitive_BQ_vs_Control_ordered,
  "Brequinar in Sensitive Cells"
)

cat("Processing Comparison 4: Interaction Effect\n")
res_interaction_annotated <- add_gene_annotations(
  res_interaction_ordered,
  "Interaction Effect"
)

# ---- Save annotated results ----

cat("SAVING ANNOTATED RESULTS\n")

write.csv(res_resistant_annotated,
          "results/DESeq2/DESeq2_Resistant_vs_Sensitive_Annotated.csv",
          row.names = FALSE)
cat("Saved: results/DESeq2/DESeq2_Resistant_vs_Sensitive_Annotated.csv\n")

write.csv(res_BQ_resistant_annotated,
          "results/DESeq2/DESeq2_Resistant_Brequinar_Annotated.csv",
          row.names = FALSE)
cat("Saved: results/DESeq2/DESeq2_Resistant_Brequinar_Annotated.csv\n")

write.csv(res_BQ_sensitive_annotated,
          "results/DESeq2/DESeq2_Sensitive_Brequinar_Annotated.csv",
          row.names = FALSE)
cat("Saved: results/DESeq2/DESeq2_Sensitive_Brequinar_Annotated.csv\n")

write.csv(res_interaction_annotated,
          "results/DESeq2/DESeq2_Interaction_Effect_Annotated.csv",
          row.names = FALSE)
cat("Saved: results/DESeq2/DESeq2_Interaction_Effect_Annotated.csv\n\n")

# ---- Display sample of annotated results ----

cat("Sample of annotated results (Resistant vs Sensitive):\n")
print(head(res_resistant_annotated[, c("gene_symbol", "ensembl_gene_id",
                                       "baseMean", "log2FoldChange", "padj", 
                                       "gene_biotype", "description")], 15))

cat("Sample of annotated results (Resistant Brequinar):\n")
print(head(res_BQ_resistant_annotated[, c("gene_symbol", "ensembl_gene_id",
                                          "baseMean", "log2FoldChange", "padj", 
                                          "gene_biotype", "description")], 15))

cat("Sample of annotated results (Sensitive Brequinar):\n")
print(head(res_BQ_sensitive_annotated[, c("gene_symbol", "ensembl_gene_id",
                                          "baseMean", "log2FoldChange", "padj", 
                                          "gene_biotype", "description")], 15))

cat("Sample of annotated results (Interaction effect):\n")
print(head(res_interaction_annotated[, c("gene_symbol", "ensembl_gene_id",
                                         "baseMean", "log2FoldChange", "padj", 
                                         "gene_biotype", "description")], 15))


# Saving annotation objects

save(res_resistant_annotated,
     res_BQ_resistant_annotated,
     res_BQ_sensitive_annotated,
     res_interaction_annotated,
     file = "results/DESeq2/DESeq2_Annotated_Results.RData")

cat("\nAll annotated results saved: results/DESeq2/DESeq2_Annotated_Results.RData\n")