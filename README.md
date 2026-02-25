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

### Results

### Discussion

## References
