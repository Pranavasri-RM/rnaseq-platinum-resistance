source("scripts/setup.R")

load("results/DESeq2/DESeq2_Results.RData")

# ============================================
# PATHWAY ENRICHMENT (GO + KEGG)
#   Up and down genes are tested separately so each pathway has a direction.
#   The background is the genes I actually tested, not the whole genome.
#   KEGG needs internet, so if it fails the GO results are still saved.
# ============================================

PADJ_CUTOFF <- 0.05
LFC_CUTOFF  <- 1
OUT_DIR     <- "results/enrichment"

cat("\nPATHWAY ENRICHMENT ANALYSIS (GO + KEGG)\n\n")

# ---- Background: all tested genes, converted to Entrez IDs ----
universe_map <- suppressMessages(
  bitr(rownames(res_resistant_vs_sensitive_ordered),
       fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
universe_entrez <- unique(universe_map$ENTREZID)
cat("Background universe:", length(universe_entrez), "genes with Entrez IDs\n\n")

to_entrez <- function(symbols) {
  if (length(symbols) == 0) return(character(0))
  unique(suppressMessages(
    bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))$ENTREZID)
}

save_enrichment <- function(obj, prefix, label, title) {
  df <- if (is.null(obj)) data.frame() else as.data.frame(obj)
  write.csv(df, file.path(OUT_DIR, paste0(prefix, "_", label, ".csv")), row.names = FALSE)
  if (nrow(df) > 0) {
    p <- dotplot(obj, showCategory = min(20, nrow(df))) + ggtitle(title)
    pdf(file.path(OUT_DIR, paste0(prefix, "_", label, "_dotplot.pdf")), width = 11, height = 9)
    print(p)
    invisible(dev.off())
  }
  nrow(df)
}

run_go <- function(genes, ont) {
  enrichGO(gene = genes, universe = universe_entrez, OrgDb = org.Hs.eg.db,
           ont = ont, pAdjustMethod = "BH",
           pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
}

run_kegg <- function(genes) {
  tryCatch({
    k <- enrichKEGG(gene = genes, universe = universe_entrez, organism = "hsa",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    if (!is.null(k)) k <- setReadable(k, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    k
  }, error = function(e) {
    cat("    KEGG unavailable:", conditionMessage(e), "\n")
    NULL
  })
}

# ---- One comparison, run separately for up and down genes ----
enrich_comparison <- function(res, comparison_name, prefix) {
  cat("===", comparison_name, "===\n")
  res <- as.data.frame(res)
  sig <- !is.na(res$padj) & res$padj < PADJ_CUTOFF

  gene_sets <- list(
    Up   = rownames(res)[sig & res$log2FoldChange >  LFC_CUTOFF],
    Down = rownames(res)[sig & res$log2FoldChange < -LFC_CUTOFF]
  )

  out <- list()
  for (dir in names(gene_sets)) {
    entrez <- to_entrez(gene_sets[[dir]])
    cat(sprintf("  %-4s: %d genes (%d with Entrez ID)\n",
                dir, length(gene_sets[[dir]]), length(entrez)))
    if (length(entrez) < 10) {
      cat("        fewer than 10 genes, skipping\n")
      out[[dir]] <- NULL
      next
    }

    lab <- paste(comparison_name, "-", dir)
    go_bp <- run_go(entrez, "BP")
    go_mf <- run_go(entrez, "MF")
    go_cc <- run_go(entrez, "CC")
    kegg  <- run_kegg(entrez)

    n <- c(
      BP   = save_enrichment(go_bp, prefix, paste0(dir, "_GO_BP"), paste(lab, ": GO Biological Process")),
      MF   = save_enrichment(go_mf, prefix, paste0(dir, "_GO_MF"), paste(lab, ": GO Molecular Function")),
      CC   = save_enrichment(go_cc, prefix, paste0(dir, "_GO_CC"), paste(lab, ": GO Cellular Component")),
      KEGG = save_enrichment(kegg,  prefix, paste0(dir, "_KEGG"),  paste(lab, ": KEGG"))
    )
    cat(sprintf("        enriched terms: GO-BP %d, GO-MF %d, GO-CC %d, KEGG %d\n",
                n["BP"], n["MF"], n["CC"], n["KEGG"]))

    out[[dir]] <- list(GO_BP = go_bp, GO_MF = go_mf, GO_CC = go_cc, KEGG = kegg,
                       genes = gene_sets[[dir]])
  }
  cat("\n")
  out
}

enrich_resistant <- enrich_comparison(res_resistant_vs_sensitive_ordered,
                                        "Resistant vs Sensitive", "Enrichment_Resistant_vs_Sensitive")
enrich_BQ_resistant <- enrich_comparison(res_resistant_BQ_vs_Control_ordered,
                                        "Brequinar in Resistant", "Enrichment_Brequinar_Resistant")
enrich_BQ_sensitive <- enrich_comparison(res_sensitive_BQ_vs_Control_ordered,
                                        "Brequinar in Sensitive", "Enrichment_Brequinar_Sensitive")
enrich_interaction <- enrich_comparison(res_interaction_ordered,
                                        "Interaction", "Enrichment_Interaction_Effect")

# ============================================
# QUICK LOOK: top GO-BP and KEGG terms (full tables are in the CSVs)
# ============================================

cat("TOP ENRICHED TERMS\n\n")

print_top <- function(enr, comparison_name, n = 5) {
  for (dir in c("Up", "Down")) {
    x <- enr[[dir]]
    if (is.null(x)) next
    for (db in c("GO_BP", "KEGG")) {
      if (is.null(x[[db]])) next
      df <- as.data.frame(x[[db]])
      if (nrow(df) == 0) next
      cat("---", comparison_name, "|", dir, "|", db, "---\n")
      print(head(df[, c("Description", "GeneRatio", "p.adjust", "Count")], n), row.names = FALSE)
      cat("\n")
    }
  }
}

print_top(enrich_resistant,    "Resistant vs Sensitive")
print_top(enrich_BQ_resistant, "Brequinar in Resistant")
print_top(enrich_BQ_sensitive, "Brequinar in Sensitive")
print_top(enrich_interaction,  "Interaction")

save(enrich_resistant, enrich_BQ_resistant, enrich_BQ_sensitive, enrich_interaction,
     universe_entrez,
     file = file.path(OUT_DIR, "Pathway_Enrichment_Results.RData"))

cat("Saved:", file.path(OUT_DIR, "Pathway_Enrichment_Results.RData"), "\n")