source("scripts/setup.R")

load("results/DESeq2/DESeq2_Results.RData")

# ============================================
# GENE ANNOTATION
#   The count files only have gene symbols, so I add Ensembl/Entrez IDs, gene
#   names and biotypes. I started with biomaRt, but Ensembl kept timing out, so
#   I switched to org.Hs.eg.db, which works offline. Symbols that map to several
#   Ensembl IDs keep one row; n_ensembl_ids says how many there were.
# ============================================

cat("\nGENE ANNOTATION\n")

all_symbols <- rownames(res_resistant_vs_sensitive_ordered)   # same 17,661 genes in every comparison

# ---- Look up all tested genes once (every comparison has the same genes) ----
cat("Annotating with org.Hs.eg.db", as.character(packageVersion("org.Hs.eg.db")), "\n")

ens_list <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = all_symbols, keytype = "SYMBOL",
                                  column = "ENSEMBL", multiVals = "list")

annotations <- data.frame(
  gene_symbol     = all_symbols,
  ensembl_gene_id = sapply(ens_list, function(x) if (all(is.na(x))) NA_character_ else sort(x)[1]),
  n_ensembl_ids   = sapply(ens_list, function(x) sum(!is.na(x))),
  entrez_id       = AnnotationDbi::mapIds(org.Hs.eg.db, all_symbols, keytype = "SYMBOL",
                                          column = "ENTREZID", multiVals = "first"),
  description     = AnnotationDbi::mapIds(org.Hs.eg.db, all_symbols, keytype = "SYMBOL",
                                          column = "GENENAME", multiVals = "first"),
  gene_biotype    = AnnotationDbi::mapIds(org.Hs.eg.db, all_symbols, keytype = "SYMBOL",
                                          column = "GENETYPE", multiVals = "first"),
  cytoband        = AnnotationDbi::mapIds(org.Hs.eg.db, all_symbols, keytype = "SYMBOL",
                                          column = "MAP", multiVals = "first"),
  row.names = NULL, stringsAsFactors = FALSE
)

cat("Symbols matched to Entrez:  ", sum(!is.na(annotations$entrez_id)), "/", length(all_symbols), "\n")
cat("Symbols matched to Ensembl: ", sum(!is.na(annotations$ensembl_gene_id)), "/", length(all_symbols), "\n")
cat("Symbols with >1 Ensembl ID: ", sum(annotations$n_ensembl_ids > 1), "\n\n")

writeLines(paste("org.Hs.eg.db version:", packageVersion("org.Hs.eg.db")),
           "results/DESeq2/annotation_version.txt")

# ---- Attach annotation to one comparison ----
annotate_results <- function(res, comparison_name) {
  df <- as.data.frame(res)
  df$gene_symbol <- rownames(df)

  df <- left_join(df, annotations, by = "gene_symbol")

  stopifnot(!anyDuplicated(df$gene_symbol), nrow(df) == nrow(res))

  # Keep genes with padj = NA instead of dropping them; they go to the bottom
  df <- df[order(df$padj, na.last = TRUE), ]
  df <- df %>% relocate(gene_symbol)

  cat(comparison_name, ":", nrow(df), "genes,",
      sum(!is.na(df$ensembl_gene_id)), "with Ensembl ID\n")
  df
}

res_resistant_annotated    <- annotate_results(res_resistant_vs_sensitive_ordered,  "Resistant vs Sensitive")
res_BQ_resistant_annotated <- annotate_results(res_resistant_BQ_vs_Control_ordered, "Brequinar in Resistant")
res_BQ_sensitive_annotated <- annotate_results(res_sensitive_BQ_vs_Control_ordered, "Brequinar in Sensitive")
res_interaction_annotated  <- annotate_results(res_interaction_ordered,             "Interaction")

# ---- Save ----
cat("\nSAVING ANNOTATED RESULTS\n")

outputs <- list(
  "results/DESeq2/DESeq2_Resistant_vs_Sensitive_Annotated.csv" = res_resistant_annotated,
  "results/DESeq2/DESeq2_Resistant_Brequinar_Annotated.csv"    = res_BQ_resistant_annotated,
  "results/DESeq2/DESeq2_Sensitive_Brequinar_Annotated.csv"    = res_BQ_sensitive_annotated,
  "results/DESeq2/DESeq2_Interaction_Effect_Annotated.csv"     = res_interaction_annotated
)
for (path in names(outputs)) {
  write.csv(outputs[[path]], path, row.names = FALSE)
  cat("Saved:", path, "\n")
}

cat("\nTop 10 annotated genes (Resistant vs Sensitive):\n")
show_cols <- intersect(c("gene_symbol", "ensembl_gene_id", "baseMean", "log2FoldChange",
                         "padj", "gene_biotype", "description"),
                       colnames(res_resistant_annotated))
print(head(res_resistant_annotated[, show_cols], 10))

save(res_resistant_annotated,
     res_BQ_resistant_annotated,
     res_BQ_sensitive_annotated,
     res_interaction_annotated,
     file = "results/DESeq2/DESeq2_Annotated_Results.RData")

cat("\nSaved: results/DESeq2/DESeq2_Annotated_Results.RData\n")