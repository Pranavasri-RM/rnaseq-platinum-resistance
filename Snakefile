rule all:
    input:
        # QC
        "results/QC/QC_PCA_Plot.pdf",

        # DESeq
        "results/DESeq2/DESeq2_Analysis_Summary.csv",

        # Visualization
        "results/visualization/Volcano_Resistant_vs_Sensitive.pdf",
        "results/visualization/MA_Plots.pdf",

        # Annotation
        "results/DESeq2/DESeq2_Annotated_Results.RData",

        # Enrichment
        "results/enrichment/Enrichment_Resistant_vs_Sensitive_GO_BP.csv",

        # Pyrimidine
        "results/enrichment/Pyrimidine_Metabolism_Summary.csv",

        # Biomarkers
        "results/biomarkers/Biomarker_Candidates_for_Brequinar_Resistance.csv"


# -------------------------------------
# 1: QC
# -------------------------------------

rule qc:
    input:
        script="scripts/01_load_and_qc.R"
    output:
        "results/QC/QC_PCA_Plot.pdf"
    shell:
        "Rscript {input.script}"


# -------------------------------------
# 2: DESEQ
# -------------------------------------

rule deseq:
    input:
        script="scripts/02_deseq_analysis.R",
        qc="results/QC/QC_PCA_Plot.pdf"
    output:
        summary="results/DESeq2/DESeq2_Analysis_Summary.csv",
        rdata="results/DESeq2/DESeq2_Results.RData"
    shell:
        "Rscript {input.script}"


# -------------------------------------
# 3: VISUALIZATION
# -------------------------------------

rule visualization:
    input:
        script="scripts/03_visualization.R",
        rdata="results/DESeq2/DESeq2_Results.RData"
    output:
        volcano="results/visualization/Volcano_Resistant_vs_Sensitive.pdf",
        ma="results/visualization/MA_Plots.pdf"
    shell:
        "Rscript {input.script}"


# -------------------------------------
# 4: ANNOTATION
# -------------------------------------

rule annotation:
    input:
        script="scripts/04_annotation.R",
        rdata="results/DESeq2/DESeq2_Results.RData"
    output:
        r_vs_s="results/DESeq2/DESeq2_Resistant_vs_Sensitive_Annotated.csv",
        bq_r="results/DESeq2/DESeq2_Resistant_Brequinar_Annotated.csv",
        bq_s="results/DESeq2/DESeq2_Sensitive_Brequinar_Annotated.csv",
        interaction="results/DESeq2/DESeq2_Interaction_Effect_Annotated.csv",
        rdata="results/DESeq2/DESeq2_Annotated_Results.RData"
    shell:
        "Rscript {input.script}"


# -------------------------------------
# 5: ENRICHMENT
# -------------------------------------

rule enrichment:
    input:
        script="scripts/05_pathway_analysis.R",
        rdata="results/DESeq2/DESeq2_Annotated_Results.RData"
    output:
        go_bp="results/enrichment/Enrichment_Resistant_vs_Sensitive_GO_BP.csv",
        rdata="results/enrichment/Pathway_Enrichment_Results.RData"
    shell:
        "Rscript {input.script}"


# -------------------------------------
# 6: PYRIMIDINE
# -------------------------------------

rule pyrimidine:
    input:
        script="scripts/06_pyrimidine_analysis.R",
        rdata="results/enrichment/Pathway_Enrichment_Results.RData"
    output:
        summary="results/enrichment/Pyrimidine_Metabolism_Summary.csv",
        rdata="results/enrichment/Pyrimidine_Metabolism_Results.RData"
    shell:
        "Rscript {input.script}"


# -------------------------------------
# 7: BIOMARKERS
# -------------------------------------

rule biomarkers:
    input:
        script="scripts/07_biomarker_analysis.R",
        rdata="results/enrichment/Pyrimidine_Metabolism_Results.RData"
    output:
        csv="results/biomarkers/Biomarker_Candidates_for_Brequinar_Resistance.csv",
        pdf="results/biomarkers/Biomarker_Candidates_Expression.pdf"
    shell:
        "Rscript {input.script}"