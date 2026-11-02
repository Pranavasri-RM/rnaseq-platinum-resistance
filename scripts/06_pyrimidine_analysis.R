source("scripts/setup.R")

load("results/DESeq2/DESeq2_Results.RData")
load("results/DESeq2/DESeq2_Annotated_Results.RData")

# ============================================
# PYRIMIDINE METABOLISM
#   Brequinar blocks DHODH, part of de novo pyrimidine synthesis. My hypothesis
#   going in was that resistant cells lean more on de novo synthesis, so here I
#   compare de novo and salvage genes directly, using a curated gene list
#   matched by exact symbol and grouped by branch.
# ============================================

PADJ_CUTOFF <- 0.05
OUT_DIR <- "results/enrichment"

pyrimidine_sets <- list(
  De_novo    = c("CAD", "DHODH", "UMPS", "CTPS1", "CTPS2", "CMPK1", "DTYMK",
                 "NME1", "NME2", "RRM1", "RRM2", "RRM2B", "TYMS", "DUT"),
  Salvage    = c("UCK1", "UCK2", "UPP1", "UPP2", "TK1", "TK2", "DCK", "DCTD",
                 "CDA", "TYMP", "NT5E", "NT5C", "NT5C2", "NT5C3A"),
  Transport  = c("SLC29A1", "SLC29A2", "SLC28A1", "SLC28A3"),
  Catabolism = c("DPYD", "DPYS", "UPB1")
)

gene_branch <- data.frame(
  gene_symbol = unlist(pyrimidine_sets, use.names = FALSE),
  branch      = rep(names(pyrimidine_sets), lengths(pyrimidine_sets)),
  stringsAsFactors = FALSE
)
gene_branch$branch <- factor(gene_branch$branch, levels = names(pyrimidine_sets))

cat("\nTARGETED ANALYSIS: PYRIMIDINE METABOLISM\n")
cat("Curated genes:", nrow(gene_branch), "across", length(pyrimidine_sets), "branches\n")

not_tested <- setdiff(gene_branch$gene_symbol, res_resistant_annotated$gene_symbol)
if (length(not_tested) > 0)
  cat("Not tested (filtered as low-count):", paste(not_tested, collapse = ", "), "\n")
cat("\n")

comparisons <- list(
  Resistant_vs_Sensitive = res_resistant_annotated,
  BQ_in_Resistant        = res_BQ_resistant_annotated,
  BQ_in_Sensitive        = res_BQ_sensitive_annotated,
  Interaction            = res_interaction_annotated
)

# ---- Gene table for each comparison ----
extract_pyrimidine_genes <- function(res_annot) {
  res_annot %>%
    inner_join(gene_branch, by = "gene_symbol") %>%
    select(branch, gene_symbol, baseMean, log2FoldChange, lfcSE, padj) %>%
    arrange(branch, padj)
}

pyrimidine_resistant    <- extract_pyrimidine_genes(res_resistant_annotated)
pyrimidine_BQ_resistant <- extract_pyrimidine_genes(res_BQ_resistant_annotated)
pyrimidine_BQ_sensitive <- extract_pyrimidine_genes(res_BQ_sensitive_annotated)
pyrimidine_interaction  <- extract_pyrimidine_genes(res_interaction_annotated)

write.csv(pyrimidine_resistant,    file.path(OUT_DIR, "Pyrimidine_Genes_Resistant_vs_Sensitive.csv"), row.names = FALSE)
write.csv(pyrimidine_BQ_resistant, file.path(OUT_DIR, "Pyrimidine_Genes_Brequinar_Resistant.csv"),   row.names = FALSE)
write.csv(pyrimidine_BQ_sensitive, file.path(OUT_DIR, "Pyrimidine_Genes_Brequinar_Sensitive.csv"),   row.names = FALSE)
write.csv(pyrimidine_interaction,  file.path(OUT_DIR, "Pyrimidine_Genes_Interaction_Effect.csv"),    row.names = FALSE)

# ---- Side-by-side table: each gene across all four comparisons ----
pick <- function(df, tag) {
  df %>% select(gene_symbol, log2FoldChange, padj) %>%
    rename_with(~ paste0(tag, "_", .x), -gene_symbol)
}
pyrimidine_wide <- gene_branch %>%
  filter(gene_symbol %in% res_resistant_annotated$gene_symbol) %>%
  left_join(pick(res_resistant_annotated,    "RvsS"),   by = "gene_symbol") %>%
  left_join(pick(res_BQ_resistant_annotated, "BQ_Res"), by = "gene_symbol") %>%
  left_join(pick(res_BQ_sensitive_annotated, "BQ_Sen"), by = "gene_symbol") %>%
  left_join(pick(res_interaction_annotated,  "Inter"),  by = "gene_symbol") %>%
  arrange(branch, gene_symbol)

write.csv(pyrimidine_wide, file.path(OUT_DIR, "Pyrimidine_Genes_All_Comparisons.csv"), row.names = FALSE)

cat("=== Resistant vs Sensitive (untreated) ===\n")
print(pyrimidine_wide %>%
        transmute(branch, gene_symbol,
                  log2FC = round(RvsS_log2FoldChange, 2),
                  padj   = signif(RvsS_padj, 2)),
      row.names = FALSE)
cat("\n")

# ---- Branch-level test ----
#   One or two significant genes are easy to over-read, so I also test each
#   branch as a set (limma::geneSetTest). It ignores correlation between genes,
#   so I treat these p-values as exploratory.
branch_summary <- function(res_annot, comparison_name) {
  stat <- sign(res_annot$log2FoldChange) * -log10(res_annot$pvalue)
  stat[is.na(stat)] <- 0

  do.call(rbind, lapply(levels(gene_branch$branch), function(b) {
    genes <- gene_branch$gene_symbol[gene_branch$branch == b]
    idx   <- res_annot$gene_symbol %in% genes
    sig   <- idx & !is.na(res_annot$padj) & res_annot$padj < PADJ_CUTOFF
    data.frame(
      Comparison    = comparison_name,
      Branch        = b,
      Genes_Tested  = sum(idx),
      Sig_Up        = sum(sig & res_annot$log2FoldChange > 0),
      Sig_Down      = sum(sig & res_annot$log2FoldChange < 0),
      Median_log2FC = round(median(res_annot$log2FoldChange[idx], na.rm = TRUE), 3),
      SetTest_p_Up   = signif(geneSetTest(idx, stat, alternative = "up"), 3),
      SetTest_p_Down = signif(geneSetTest(idx, stat, alternative = "down"), 3)
    )
  }))
}

summary_pyrimidine <- do.call(rbind, Map(branch_summary, comparisons, names(comparisons)))
rownames(summary_pyrimidine) <- NULL

cat("=== Branch summary (padj < 0.05 counts; rank-based set test) ===\n")
print(summary_pyrimidine, row.names = FALSE)
cat("\n")

write.csv(summary_pyrimidine, file.path(OUT_DIR, "Pyrimidine_Metabolism_Summary.csv"), row.names = FALSE)

# ---- DHODH, brequinar's direct target ----
cat("=== DHODH across comparisons ===\n")
print(pyrimidine_wide %>% filter(gene_symbol == "DHODH") %>% select(-branch), row.names = FALSE)
cat("\n")

# ---- Which genes are behind the pyrimidine GO term from step 05? ----
go_file <- file.path(OUT_DIR, "Enrichment_Resistant_vs_Sensitive_Down_GO_BP.csv")
if (file.exists(go_file)) {
  go <- read.csv(go_file)
  hit <- go[grepl("pyrimidine", go$Description, ignore.case = TRUE), c("Description", "p.adjust", "geneID")]
  if (nrow(hit) > 0) {
    cat("=== Pyrimidine-related GO term from step 05 (down in resistant) ===\n")
    print(hit, row.names = FALSE)
    cat("\n")
  }
}

# ---- Heatmap of all tested pyrimidine genes (VST, row z-scores), grouped by branch ----
vsd <- vst(dds, blind = FALSE)
genes_hm <- pyrimidine_wide$gene_symbol
z <- t(scale(t(assay(vsd)[genes_hm, , drop = FALSE])))
z[z > 3] <- 3; z[z < -3] <- -3

row_anno <- data.frame(Branch = pyrimidine_wide$branch, row.names = genes_hm)
col_anno <- data.frame(Resistance = sample_info$Resistance,
                       Treatment  = sample_info$Treatment,
                       row.names  = rownames(sample_info))
anno_colors <- list(
  Resistance = c(Sensitive = "blue", Resistant = "red"),
  Treatment  = c(Control = "gray70", Brequinar = "orange"),
  Branch     = c(De_novo = "#1b9e77", Salvage = "#d95f02", Transport = "#7570b3", Catabolism = "#e7298a")
)

pdf(file.path(OUT_DIR, "Heatmap_Pyrimidine_Genes.pdf"), width = 9, height = 10)
pheatmap(z,
         annotation_row = row_anno, annotation_col = col_anno,
         annotation_colors = anno_colors,
         cluster_rows = FALSE, cluster_cols = TRUE,
         gaps_row = cumsum(table(pyrimidine_wide$branch))[-length(pyrimidine_sets)],
         main = "Pyrimidine metabolism genes (VST, row z-score)",
         fontsize_row = 8,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-3, 3, length.out = 101))
invisible(dev.off())
cat("Saved:", file.path(OUT_DIR, "Heatmap_Pyrimidine_Genes.pdf"), "\n")

# ---- Save ----
save(pyrimidine_resistant, pyrimidine_BQ_resistant, pyrimidine_BQ_sensitive,
     pyrimidine_interaction, pyrimidine_wide, summary_pyrimidine, pyrimidine_sets,
     file = file.path(OUT_DIR, "Pyrimidine_Metabolism_Results.RData"))

cat("Saved:", file.path(OUT_DIR, "Pyrimidine_Metabolism_Results.RData"), "\n")