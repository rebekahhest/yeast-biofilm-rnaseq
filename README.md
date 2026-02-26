# yeast-biofilm-rnaseq
Differential gene expression analysis and functional enrichment/annotation analysis for biofilm maturation in flor yeast, Saccharomyces cerevisiae

**Author:** Rebekah Hest  
**Course:** BINF6110 - Genomic Methods for Bioinformatics

## Assignment 2

This repository contains Assignment 2, a study to examine differential expression in a bulk transcriptomics experiment and perform functional annotations for yeast biofilm (velum) development during wine aging. 

### Introduction

Biofilm formation is a part of a complex, dynamic aging process during sherry wine production (Carbonero-Pacheco et al., 2022). A strain of flor yeast, _Saccharomyces cerevisiae_, is responsible for the biochemical changes needed to survive the stressful conditions of biological aging (Mardanov et al., 2020). The structure formed on the air-liquid interface of the wine (flor velum) acts as a protective barrier in oxidative metabolism and contributes to the aroma and flavour development (Moreno-Garcia et al., 2017). The developmental transition from early biofilm to mature velum reflects a robust, specialized response to environmental and nutritional changes.
<br><br>
	The adaptations underpinning the regulation of flor formation can be characterized through the expression of specific genes (Legras et al., 2016).  Bulk RNA sequencing (RNA-seq) is a comprehensive and cost-effective approach for quantifying transcript abundance and identifying genes that are differentially expressed between experimental conditions (Wang et al., 2009). Comparison of transcriptional profiles across early, thin, and mature biofilm stages can reveal the molecular pathways and regulatory networks underlying velum development. However, detecting transcriptome changes requires a rigorous statistical framework for count-based RNA-seq data to account for biological variability, technical noise, and handle multiple hypothesis testing (Love et al., 2014).
  <br><br>
Several computational tools have been developed for differential expression analysis (DEA), including DESeq2, edgeR, and limma-voom (Kalantari-Dehaghi et al, 2025). DESeq2 and edgeR model count data using the negative binomial distribution, whereas limma-voom estimates precision weights within a linear modeling framework (Robinson et al., 2009; Law et al., 2014). However, DESeq2 incorporates shrinkage estimations for dispersion and fold changes to improve statistical stability in experiments with small to moderate sample sizes (Love et al., 2014).
<br><br>
Pathway analysis can be performed using overrepresentation analysis (ORA) or gene set enrichment analysis (GSEA). ORA evaluates whether predefined gene sets are significantly expressed between conditions in a hypergeometric distribution and assumes genes in a pathway are independent of each other (Ashburner et al., 2000; Khatri et al., 2012). In contrast, GSEA analyzes a ranked gene list to detect subtle, pathway-level changes (Subramanian et al., 2005). Because this dataset includes clear stage-specific contrasts with robust differential expression, ORA was selected to identify significantly enriched biological processes in a straightforward and interpretable manner.

### Methods

### Data Acquisition

RNA-seq data were obtained from the NCBI BioProject PRJNA592304, which includes transcriptomic profiles of _Saccharomyces cerevisiae_ during early, thin, and mature biofilm (velum) development. Three biological replicates were available per developmental stage. Raw sequencing reads were downloaded from the Sequence Read Archive (SRA) using the SRA Toolkit (v3.2.1) (`prefetch` and `fasterq-dump`). The _Saccharomyces cerevisiae_ reference transcriptome (GCF_000146045.2_R64) was retrieved from the NCBI RefSeq database in FASTA format for transcript quantification. 

### Quality Control

Quality assessment of raw sequences was performed using FastQC (v0.12.1) to evaluate read quality scores, lengths and GC content (`fastqc`). No trimming was performed to reduce the noted impact of aggressive base removal on gene expression estimates (Williams et al., 2016).
	
### Transcriptome Quantification

RNA-seq quantification was performed using Salmon (v1.10.34), a pseudoalignment tool that accurately estimates transcript abundance and models biases comparable to splice-aware alignment methods while substantial reducing computational time and storage requirements (Patro et al., 2017; Schaarschmidt et al., 2020). A transcriptome index was constructed from the reference (`index`) with default k-mer length parameters. Transcripts were quantified independently using quasi-mapping (`salmon quant`) command with automatic library selection (`-l A`) and selective alignment (`--validateMappings`) against the index library.
	
### Differential Expression Analysis

Transcript-level abundances estimated using Salmon were imported into R (v2026.01.1+403) using the `tximport` package (v1.36.1). Quantification data was converted to gene-level counts using a mapping derived from the reference GTF annotation. Differential expression analysis (DEA) was performed using DESeq2 (v1.48.2) to fit a negative binomial GLM with biofilm stage as the condition (`DEQeq`). Wald tests were applied to estimate log2 fold changes between conditions, and p-values were adjusted for multiple testing using the Benjamini–Hochberg false discovery rate method (Love et al., 2014). Abundance estimates were normalization via Variance stabilizing transformation (VST) (`vst`). Principal component analysis (PCA) was performed to assess global differences between biofilm stages. Shrinkage was applied to the normalized DEA results using the apeglm method (`lfcShrink`) and MA plots were generated to display log fold change relative to mean expression counts. Volcano plots were generated to visualize significantly up- and downregulated genes (adjusted p-value < 0.05; log2 fold change > 1). A heatmap was produced for the top 20 most significantly expressed genes using the pheatmap package (v1.0.13).

### Functional Annotation

Gene Ontology (GO) enrichment analysis was performed using the clusterProfiler package (v4.16.0) for overrepresentation analysis (ORA) to identify biological processes (BP) (`enrichGO`). The Saccharomyces Genome Database’s (SGD’s) was used for annotating (org.Sc.sgd.db) where the open reading frames (ORF) (or systematic names) were used as the gene IDs (Dwight et al., 2002). The background gene universe was composed of all unique genes modelled in the DEA. Statistical significance was calculated using the hypergeometric distribution with Benjamini–Hochberg (BH) correction (adjusted p-value < 0.05). Enriched GO terms were visualized using a dot plot. 

### Results

### Discussion

## References
