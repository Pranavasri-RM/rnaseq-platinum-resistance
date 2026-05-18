source("scripts/setup.R")

# From DESeq step
load("results/DESeq2/DESeq2_Results.RData")

# From Annotation step
load("results/DESeq2/DESeq2_Annotated_Results.RData")

# ============================================
#  PATHWAY ENRICHMENT ANALYSIS
# ============================================

cat("\nPATHWAY ENRICHMENT ANALYSIS\n")
cat("Gene Ontology (GO) & KEGG Pathways\n")

# Function for GO and KEGG enrichment
perform_pathway_enrichment <- function(res, comparison_name, title_prefix) {
  
  cat("Processing:", comparison_name, "\n")
  
  # Get significant genes (padj < 0.05 AND |log2FC| > 1)
  sig_indices <- which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)
  sig_genes <- rownames(res)[sig_indices]
  
  n_sig <- length(sig_genes)
  cat("Found", n_sig, "significant DEGs\n")
  
  if (n_sig == 0) {
    cat("No significant genes found for enrichment\n\n")
    return(NULL)
  }
  
  tryCatch({
    # Convert gene symbols to Entrez IDs
    cat("  Converting to Entrez IDs...\n")
    gene_list <- bitr(sig_genes,
                      fromType = "SYMBOL",        # Changed from ENSEMBL to SYMBOL
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)
    
    n_converted <- nrow(gene_list)
    cat("  Successfully converted:", n_converted, "/", n_sig, "genes\n")
    
    if (n_converted == 0) {
      cat("No genes could be converted to Entrez IDs\n\n")
      return(NULL)
    }
    
    # ---- GO ENRICHMENT - BIOLOGICAL PROCESS ----
    cat("Performing GO Biological Process enrichment...\n")
    ego_BP <- enrichGO(gene = gene_list$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
    
    n_BP <- nrow(as.data.frame(ego_BP))
    cat("    Found", n_BP, "enriched BP terms\n")
    
    # ---- GO ENRICHMENT - MOLECULAR FUNCTION ----
    cat("Performing GO Molecular Function enrichment...\n")
    ego_MF <- enrichGO(gene = gene_list$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "MF",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
    
    n_MF <- nrow(as.data.frame(ego_MF))
    cat("Found", n_MF, "enriched MF terms\n")
    
    # ---- GO ENRICHMENT - CELLULAR COMPONENT ----
    cat("  Performing GO Cellular Component enrichment...\n")
    ego_CC <- enrichGO(gene = gene_list$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05,
                       readable = TRUE)
    
    n_CC <- nrow(as.data.frame(ego_CC))
    cat("    Found", n_CC, "enriched CC terms\n")
    
    # ---- KEGG ENRICHMENT ----
    cat("  Performing KEGG pathway enrichment...\n")
    kegg <- enrichKEGG(gene = gene_list$ENTREZID,
                       organism = "hsa",
                       pvalueCutoff = 0.05)
    
    n_KEGG <- nrow(as.data.frame(kegg))
    cat("    Found", n_KEGG, "enriched KEGG pathways\n\n")
    
    # ---- VISUALIZATIONS ----
    
    # GO Biological Process
    if (n_BP > 0) {
      cat("  Creating GO BP visualizations...\n")
      
      # Dot plot
      pdf(paste0("results/enrichment/", title_prefix, "_GO_BP_dotplot.pdf"), width = 12, height = 10)
      p1 <- dotplot(ego_BP, showCategory = min(20, n_BP)) + 
        ggtitle(paste(comparison_name, "- GO Biological Process"))
      print(p1)
      dev.off()
      cat("Saved:", paste0("results/enrichment/", title_prefix, "_GO_BP_dotplot.pdf\n"))
      
      # Bar plot
      pdf(paste0("results/enrichment/", title_prefix, "_GO_BP_barplot.pdf"), width = 12, height = 10)
      p2 <- barplot(ego_BP, showCategory = min(20, n_BP)) +
        ggtitle(paste(comparison_name, "- GO Biological Process"))
      print(p2)
      dev.off()
      cat("Saved:", paste0("results/enrichment/", title_prefix, "_GO_BP_barplot.pdf\n"))
    } else {
      cat("No BP terms to visualize\n")
    }
    
    # GO Molecular Function
    if (n_MF > 0) {
      cat("Creating GO MF visualizations...\n")
      
      pdf(paste0("results/enrichment/", title_prefix, "_GO_MF_dotplot.pdf"), width = 12, height = 10)
      p3 <- dotplot(ego_MF, showCategory = min(20, n_MF)) +
        ggtitle(paste(comparison_name, "- GO Molecular Function"))
      print(p3)
      dev.off()
      cat("Saved:", paste0("results/enrichment/", title_prefix, "_GO_MF_dotplot.pdf\n"))
    } else {
      cat("No MF terms to visualize\n")
    }
    
    # KEGG Pathways
    if (n_KEGG > 0) {
      cat("Creating KEGG visualizations...\n")
      
      pdf(paste0("results/enrichment/", title_prefix, "_KEGG_dotplot.pdf"), width = 12, height = 10)
      p4 <- dotplot(kegg, showCategory = min(20, n_KEGG)) +
        ggtitle(paste(comparison_name, "- KEGG Pathways"))
      print(p4)
      dev.off()
      cat("Saved:", paste0("results/enrichment/", title_prefix, "_KEGG_dotplot.pdf\n"))
    } else {
      cat("No KEGG pathways to visualize\n")
    }
    
    # ---- SAVE RESULTS ----
    cat("Saving enrichment results...\n")
    
    write.csv(as.data.frame(ego_BP), 
              paste0("results/enrichment/", title_prefix, "_GO_BP.csv"), row.names = FALSE)
    cat("Saved:", paste0("results/enrichment/", title_prefix, "_GO_BP.csv\n"))
    
    write.csv(as.data.frame(ego_MF),
              paste0("results/enrichment/", title_prefix, "_GO_MF.csv"), row.names = FALSE)
    cat("Saved:", paste0("results/enrichment/", title_prefix, "_GO_MF.csv\n"))
    
    write.csv(as.data.frame(ego_CC),
              paste0("results/enrichment/", title_prefix, "_GO_CC.csv"), row.names = FALSE)
    cat("Saved:", paste0("results/enrichment/", title_prefix, "_GO_CC.csv\n"))
    
    write.csv(as.data.frame(kegg),
              paste0("results/enrichment/", title_prefix, "_KEGG.csv"), row.names = FALSE)
    cat("Saved:", paste0("results/enrichment/", title_prefix, "_KEGG.csv\n\n"))
    
    # Return all results
    return(list(GO_BP = ego_BP, 
                GO_MF = ego_MF,
                GO_CC = ego_CC,
                KEGG = kegg,
                gene_list = gene_list))
    
  }, error = function(e) {
    cat("Error during enrichment:", conditionMessage(e), "\n\n")
    return(NULL)
  })
}

# ============================================
# RUN ENRICHMENT FOR ALL COMPARISONS
# ============================================

cat("Running pathway enrichment for all comparisons...\n\n")

enrich_resistant <- perform_pathway_enrichment(
  res_resistant_vs_sensitive_ordered,
  "Resistant vs Sensitive (Control)",
  "Enrichment_Resistant_vs_Sensitive"
)

enrich_BQ_resistant <- perform_pathway_enrichment(
  res_resistant_BQ_vs_Control_ordered,
  "Brequinar in Resistant Cells",
  "Enrichment_Brequinar_Resistant"
)

enrich_BQ_sensitive <- perform_pathway_enrichment(
  res_sensitive_BQ_vs_Control_ordered,
  "Brequinar in Sensitive Cells",
  "Enrichment_Brequinar_Sensitive"
)

enrich_interaction <- perform_pathway_enrichment(
  res_interaction_ordered,
  "Interaction Effect",
  "Enrichment_Interaction_Effect"
)

# ============================================
# SUMMARY AND TOP PATHWAYS
# ============================================

cat("ENRICHMENT ANALYSIS SUMMARY\n")

print_top_pathways <- function(enrichment_list, comparison_name) {
  if (is.null(enrichment_list)) {
    cat(comparison_name, "- No significant pathways\n\n")
    return()
  }
  
  cat("---", comparison_name, "---\n")
  
  if (!is.null(enrichment_list$GO_BP)) {
    df_BP <- as.data.frame(enrichment_list$GO_BP)
    if (nrow(df_BP) > 0) {
      cat("Top 5 GO Biological Processes:\n")
      top_BP <- head(df_BP[, c("ID", "Description", "p.adjust", "Count")], 5)
      print(top_BP)
      cat("\n")
    }
  }
  
  if (!is.null(enrichment_list$KEGG)) {
    df_KEGG <- as.data.frame(enrichment_list$KEGG)
    if (nrow(df_KEGG) > 0) {
      cat("Top 5 KEGG Pathways:\n")
      top_KEGG <- head(df_KEGG[, c("ID", "Description", "p.adjust", "Count")], 5)
      print(top_KEGG)
      cat("\n")
    }
  }
  
  cat("\n")
}

print_top_pathways(enrich_resistant, "Resistant vs Sensitive")
print_top_pathways(enrich_BQ_resistant, "Brequinar in Resistant")
print_top_pathways(enrich_BQ_sensitive, "Brequinar in Sensitive")
print_top_pathways(enrich_interaction, "Interaction Effect")

# ============================================
# SAVE ALL ENRICHMENT OBJECTS
# ============================================

save(enrich_resistant, 
     enrich_BQ_resistant, 
     enrich_BQ_sensitive,
     enrich_interaction,
     file = "results/enrichment/Pathway_Enrichment_Results.RData")

cat("All enrichment results saved: results/enrichment/Pathway_Enrichment_Results.RData\n")