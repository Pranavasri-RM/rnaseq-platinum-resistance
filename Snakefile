import glob

# Raw GEO count files (GSM8661174-GSM8661185)
COUNTS = sorted(glob.glob("data/*.count.out.txt"))
SETUP  = "scripts/setup.R"

rule all:
    input:
        "results/QC/QC_PCA_Plot.pdf",
        "results/DESeq2/DESeq2_Analysis_Summary.csv",
        "results/visualization/Volcano_Resistant_vs_Sensitive.pdf",
        "results/DESeq2/DESeq2_Annotated_Results.RData",
        "results/enrichment/Pathway_Enrichment_Results.RData",
        "results/enrichment/Pyrimidine_Metabolism_Summary.csv",
        "results/biomarkers/Resistance_Marker_Candidates.csv",
        "results/biomarkers/Final_Pipeline_Summary.csv"


# 1: Load counts, filter, TMM log-CPM QC (PCA, distance and correlation heatmaps)
rule qc:
    input:
        script="scripts/01_load_and_qc.R",
        setup=SETUP,
        counts=COUNTS
    output:
        pca="results/QC/QC_PCA_Plot.pdf",
        pca_png="results/figures/pca_plot.png",
        rdata="results/QC/Filtered_Data.RData"
    shell:
        "Rscript {input.script}"


# 2: DESeq2 (3 pairwise contrasts + interaction model), ashr-shrunken fold changes
rule deseq:
    input:
        script="scripts/02_deseq_analysis.R",
        setup=SETUP,
        rdata="results/QC/Filtered_Data.RData"
    output:
        summary="results/DESeq2/DESeq2_Analysis_Summary.csv",
        rdata="results/DESeq2/DESeq2_Results.RData"
    shell:
        "Rscript {input.script}"


# 3: Volcano, MA and top-50 heatmaps (VST)
rule visualization:
    input:
        script="scripts/03_visualization.R",
        setup=SETUP,
        rdata="results/DESeq2/DESeq2_Results.RData"
    output:
        volcano="results/visualization/Volcano_Resistant_vs_Sensitive.pdf",
        volcano_png="results/figures/volcano_resistant_vs_sensitive.png",
        ma="results/visualization/MA_Plots.pdf"
    shell:
        "Rscript {input.script}"


# 4: Gene annotation (org.Hs.eg.db, offline)
rule annotation:
    input:
        script="scripts/04_annotation.R",
        setup=SETUP,
        rdata="results/DESeq2/DESeq2_Results.RData"
    output:
        r_vs_s="results/DESeq2/DESeq2_Resistant_vs_Sensitive_Annotated.csv",
        rdata="results/DESeq2/DESeq2_Annotated_Results.RData"
    shell:
        "Rscript {input.script}"


# 5: GO + KEGG over-representation, up and down genes separately (KEGG needs internet)
rule enrichment:
    input:
        script="scripts/05_pathway_analysis.R",
        setup=SETUP,
        rdata="results/DESeq2/DESeq2_Results.RData"
    output:
        up="results/enrichment/Enrichment_Resistant_vs_Sensitive_Up_GO_BP.csv",
        down="results/enrichment/Enrichment_Resistant_vs_Sensitive_Down_GO_BP.csv",
        rdata="results/enrichment/Pathway_Enrichment_Results.RData"
    shell:
        "Rscript {input.script}"


# 6: Targeted pyrimidine metabolism analysis (curated gene list, branch-level tests)
rule pyrimidine:
    input:
        script="scripts/06_pyrimidine_analysis.R",
        setup=SETUP,
        deseq="results/DESeq2/DESeq2_Results.RData",
        annot="results/DESeq2/DESeq2_Annotated_Results.RData",
        go_down="results/enrichment/Enrichment_Resistant_vs_Sensitive_Down_GO_BP.csv"
    output:
        summary="results/enrichment/Pyrimidine_Metabolism_Summary.csv",
        genes="results/enrichment/Pyrimidine_Genes_All_Comparisons.csv",
        heatmap="results/enrichment/Heatmap_Pyrimidine_Genes.pdf",
        rdata="results/enrichment/Pyrimidine_Metabolism_Results.RData"
    shell:
        "Rscript {input.script}"


# 7: Candidate resistance markers + final summary + sessionInfo
rule markers:
    input:
        script="scripts/07_biomarker_analysis.R",
        setup=SETUP,
        deseq="results/DESeq2/DESeq2_Results.RData",
        annot="results/DESeq2/DESeq2_Annotated_Results.RData",
        pyr="results/enrichment/Pyrimidine_Metabolism_Results.RData"
    output:
        csv="results/biomarkers/Resistance_Marker_Candidates.csv",
        pdf="results/biomarkers/Resistance_Marker_Candidates_Expression.pdf",
        summary="results/biomarkers/Final_Pipeline_Summary.csv",
        session="results/sessionInfo.txt"
    shell:
        "Rscript {input.script}"
