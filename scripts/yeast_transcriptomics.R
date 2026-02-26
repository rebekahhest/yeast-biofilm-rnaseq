#### Load Libraries ------
BiocManager::install(c("tximport", "DESeq2", "readr","clusterProfiler","org.Sc.sgd.db"))
library(tximport)
library(DESeq2)
library(readr)
library(GenomicFeatures)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(clusterProfiler)
library(org.Sc.sgd.db)

#### Prepare Yeast Dataset ------

# Create Sample Table
sampleTable <- data.frame(
  Run = c("SRR10551665","SRR10551664","SRR10551663",
          "SRR10551662","SRR10551661","SRR10551660",
          "SRR10551659","SRR10551658","SRR10551657"),
  condition = c(rep("early",3),
                rep("thin",3),
                rep("mature",3))
)

# Ensure row names are exact
rownames(sampleTable) <- sampleTable$Run

# Create condition is a factor of development stages
sampleTable$condition <- factor(sampleTable$condition, levels = c("early","thin","mature"))

# Create a Transcript Database (TxDb) object from the GTF file
txdb <- makeTxDbFromGFF("data/reference/GCF_000146045.2_R64_genomic.gtf")

# Get all transcript names/IDs
k <- keys(txdb, keytype = "TXNAME")

# Select the transcript ID and gene ID columns
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# Create files vector
files <- file.path("data/quant", sampleTable$Run, "quant.sf")

# Import transcript-level abundances
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

#### Differential Gene Expression ------

# DESeq2
dds <- DESeqDataSetFromTximport(txi,colData = sampleTable, design = ~ condition)
dds <- DESeq(dds)

# Review comparisons
resultsNames(dds)

# Retrive results
res_thin_early <- results(dds, name="condition_thin_vs_early")
res_mature_early <- results(dds, name="condition_mature_vs_early")
res_mature_thin <- results(dds, contrast = c("condition", "mature", "thin"))

# Compare number of significantly expressed genes (FDR < 0.05)
sum(res_thin_early$padj < 0.05, na.rm=TRUE)
# 2202
sum(res_mature_early$padj < 0.05, na.rm=TRUE)
# 2993
sum(res_mature_thin$padj < 0.05, na.rm=TRUE)
#2374

# Shrinkage

# Need to relevel the reference as thin for shrinkage since apeglm only uses 'coef' (cannot use contrast)
dds_thin_ref <- dds
dds_thin_ref$condition <- relevel(dds_thin_ref$condition, ref="thin")
dds_thin_ref <- DESeq(dds_thin_ref)

resLFC_thin_early <- lfcShrink(dds, coef="condition_thin_vs_early", type="apeglm")
resLFC_mature_early <- lfcShrink(dds, coef="condition_mature_vs_early", type="apeglm")
resLFC_mature_thin <- lfcShrink(dds_thin_ref, coef="condition_mature_vs_thin", type="apeglm")

# MA plots with shrinkage
plotMA(resLFC_thin_early, ylim=c(-2,2))
plotMA(resLFC_mature_early, ylim=c(-2,2))
plotMA(resLFC_mature_thin, ylim=c(-2,2))

# Volcano Plot
# mark genes as upregulated, downregulated, or not significant
# create a function exclude genes with under 2-fold change (log2FoldChange < 1) for each condition comparison
make_volcano_df <- function(res_obj, contrast_name) {
  df <- as.data.frame(res_obj)
  df$gene <- rownames(df)
  
  df$significant <- ifelse(df$padj < 0.05 & abs(df$log2FoldChange) > 1,
                           ifelse(df$log2FoldChange > 0, "Up", "Down"),
                           "Not Sig")
  
  df$contrast <- contrast_name
  
  # Remove NA values
  df <- na.omit(df)
  return(df)
}

df1 <- make_volcano_df(resLFC_thin_early, "Early vs Thin")
df2 <- make_volcano_df(resLFC_mature_thin, "Thin vs Mature")
df3 <- make_volcano_df(resLFC_mature_early, "Early vs Mature")

# Retrieve count of upregulated, downregulated and non-signifcant genes
table(df1$significant)
# 485 upregulated; 453 downregulated
table(df2$significant)
# 625 upregulated; 466 downregulated 
table(df3$significant)
# 947 upregulated; 769 downregulated

# Combine the data frames to plot 3 condition comparisons together
volcano_df <- bind_rows(df1, df2, df3)

# Force condition order so comparisons reflect progressive development
volcano_df$contrast <- factor(
  volcano_df$contrast,
  levels = c("Early vs Thin",
             "Thin vs Mature",
             "Early vs Mature")
)

ggplot(volcano_df,aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point(alpha=0.5, size=1) +
  scale_color_manual(values=c("Down"="blue", "Not Sig"="gray", "Up"="red")) +
  facet_wrap(~contrast) +
  labs(x="Log2 Fold Change", y="-Log10 Adjusted p-value") +
  theme_minimal()

# Select top 20 genes by p values from non-NA genes for mature versus early pairwise comparison
top_genes <- head(order(resLFC_mature_early$padj), 20)
gene_names <- rownames(resLFC_mature_early)[top_genes]

# Normalize counts via a variance stabilizing transformation
vsd <- vst(dds)
# Store counts in a matrix for the heatmap
mat <- assay(vsd)[gene_names, ]

# Add annotation for biofilm stage (condition)
annotation_df <- data.frame(condition=sampleTable$condition)
rownames(annotation_df) <- sampleTable$Run

# Heatmap
pheatmap(mat,
         scale="row",
         cluster_rows=TRUE,
         cluster_cols=TRUE,
         annotation_col=annotation_df,
         show_rownames=TRUE,
         show_colnames=FALSE)

# Same VST used for the heatmap
vsd <- vst(dds)

# Retrieve the coordinates using plotPCA from DESeq2
pca_data <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
summary(pca_data)

# Calculate percent variance explained by the top two principal components
percentVar <- round(100 * attr(pca_data,"percentVar"))

# Plot stage by colour 
ggplot(pca_data,
       aes(PC1, PC2, color=condition, shape = condition)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of Biofilm Stages") +
  coord_fixed()

#### Functional Enrichment Analysis ------

# create a function to run the enrichment analysis for each pairwise comparison
run_go_enrichment <- function(res_object, contrast_name) {
  
  res_df <- as.data.frame(res_object)
  res_df$ORF <- rownames(res_df)
  
  # Significant genes
  sig_genes <- res_df %>%
    filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
    pull(ORF) %>%
    na.omit() %>%
    unique()
  
  # Background universe
  all_genes <- res_df %>%
    pull(ORF) %>%
    na.omit() %>%
    unique()
  
  # Enrichment
  ego <- enrichGO(
    gene          = sig_genes,
    universe      = all_genes,
    OrgDb         = org.Sc.sgd.db,
    keyType       = "ORF",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05
  )
  
  ego_df <- as.data.frame(ego)
  ego_df$contrast <- contrast_name
  
  return(list(ego_object = ego, ego_df = ego_df))
}

go_thin_early   <- run_go_enrichment(res_thin_early, "Thin vs Early")
go_mature_thin  <- run_go_enrichment(res_mature_thin, "Mature vs Thin")
go_mature_early <- run_go_enrichment(res_mature_early, "Mature vs Early")

p1 <- dotplot(go_thin_early$ego_object, showCategory=10) + ggtitle("Thin vs Early")
p2 <- dotplot(go_mature_thin$ego_object, showCategory=10) + ggtitle("Mature vs Thin")
p3 <- dotplot(go_mature_early$ego_object, showCategory=10) + ggtitle("Mature vs Early")

library(patchwork)
p1 + p2 + p3

combined_plot <- p1 + p2 + p3

ggsave(
  "GO_enrichment_combined.png",
  plot = combined_plot,
  width = 16, 
  height = 6, 
  dpi = 300
)
