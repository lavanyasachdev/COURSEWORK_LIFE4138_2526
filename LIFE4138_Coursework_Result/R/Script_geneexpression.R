#reading in the data
df1 <- read.delim("/Users/apple/Desktop/Bioinformatics Modules/Assessment/Computational coursework/Coursework/Datasets/A_vs_C.deseq2.results.tsv")
head(df1)

df2 <- read_delim('/Users/apple/Desktop/Bioinformatics Modules/Assessment/Computational coursework/Coursework/Datasets/A_vs_D.deseq2.results.tsv')
head(df2)

#loading package
library(dbplyr)
library(tidyverse)

#removing NA padj values for df1 i.e. A_vs_C dataset
df1_filtered <- df1 %>%
  filter(!is.na(padj))
print(paste("filtered genes:", nrow(df1_filtered)))

#finding significant genes for A_vs_C 
df1_sig <- df1_filtered %>%
  filter(padj < 0.05,
         abs(log2FoldChange) > 1)
print(paste("significant genes:", nrow(df1_sig)))

#finding upregulated genes for A_vs_C 
df1_up <- df1_sig %>%
  filter(log2FoldChange > 0)
print(paste("upregulated genes:", nrow(df1_up)))

#finding downregulated genes for A_vs_C 
df1_down <- df1_sig %>%
  filter(log2FoldChange < 0)
print(paste("downregulated genes:", nrow(df1_down)))

#summary of pvalues for A_vs_C
df1_p_summary <- summary(df1_filtered$pvalue)
print("summary of pvalues:") 
print(df1_p_summary)

#summary of log fold changes for A_vs_C
df1_log_summary <- summary(df1_filtered$log2FoldChange)
print("summary of log2Foldchange:") 
print(df1_log_summary)

#Dataset A_vs_D

#removing NA padj values for df2 i.e. A_vs_D dataset
df2_filtered <- df2 %>%
  filter(!is.na(padj))
print(paste("filtered genes:", nrow(df2_filtered)))

#finding significant genes for A_vs_D 
df2_sig <- df2_filtered %>%
  filter(padj < 0.05,
         abs(log2FoldChange) > 1)
print(paste("significant genes:", nrow(df2_sig)))

#finding upregulated genes for A_vs_D
df2_up <- df2_sig %>%
  filter(log2FoldChange > 0)
print(paste("upregulated genes:", nrow(df2_up)))

#finding downregulated genes for A_vs_D
df2_down <- df2_sig %>%
  filter(log2FoldChange < 0)
print(paste("downregulated genes:", nrow(df2_down)))

#summary of pvalues for A_vs_D
df2_p_summary <- summary(df2_filtered$pvalue)
print("summary of pvalues:") 
print(df2_p_summary)

#summary of log fold changes for A_vs_D
df2_log_summary <- summary(df2_filtered$log2FoldChange)
print("summary of log2Foldchange:") 
print(df2_log_summary)

#Plots

#1. Volcano plot

library(tidyverse)
library(ggplot2)

install.packages("ggrepel")
library(ggrepel)

#for A vs C
df1 <- df1 %>%
  mutate(status = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
    padj < 0.05 & log2FoldChange < 0 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

ggplot(df1, aes(log2FoldChange, -log10(pvalue), color = status)) +
  geom_point() +
  theme_light() +
  labs(title = "Volcano plot of A vs C")

# for A vs D
df2 <- df2 %>%
  mutate(status = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
    padj < 0.05 & log2FoldChange < 0 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

ggplot(df2, aes(log2FoldChange, -log10(pvalue), color = status)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Volcano plot of A vs D")

# MA plot for A vs C 

install.packages("ggpubr")
library(ggpubr)
ggmaplot(df1, fdr = 0.05, fc = 1, genenames = df1$gene_id, main = "MA plot: A vs C")

#MA plot for A vs D 
ggmaplot(df2, fdr = 0.05, fc = 1, genenames = df2$gene_id, main = "MA plot: A vs D")

#Histogram A vs C

ggplot(df1, aes(pvalue)) +
  geom_histogram (bins = 50) +
  theme_minimal() +
  labs(
    title = "Histogram of p-values of A vs C",
    x = "p-value", 
    y = "count"
  )

#Histogram A vs D

ggplot(df2, aes(pvalue)) +
  geom_histogram (bins = 50) +
  theme_minimal() +
  labs(
    title = "Histogram of p-values of A vs D",
    x = "p-value", 
    y = "count"
  )
         
# Heatmap A vs C

install.packages("pheatmap")
library(pheatmap)

# selecting top 20 differentially expressed genes
top20_1 <- df1 %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice(1:20)
#creating matrix for pheatmap
matrix1 <- as.matrix(top20_1$log2FoldChange)
rownames(matrix1) <- top20_1$gene_id

pheatmap(
  matrix1, 
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Top 20 differentially epresssed genes for A vs C",
  color = colorRampPalette(c("blue", "yellow", "pink"))(100)
)

# Heatmap A vs D

# selecting top 20 differentially expressed genes
top20_2 <- df2 %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice(1:20)
#creating matrix for pheatmap
matrix2 <- as.matrix(top20_1$log2FoldChange)
rownames(matrix2) <- top20_1$gene_id

pheatmap(
  matrix2, 
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Top 20 differentially expresssed genes for A vs D",
  color = colorRampPalette(c("red", "blue", "white"))(100)
)

#making a combined heatmap for both datasets

#finding common genes
common_genes <- intersect(df1$gene_id, df2$gene_id)

df1_common <- df1 %>% filter(gene_id %in% common_genes)
df2_common <- df2 %>% filter(gene_id %in% common_genes)

top20_common <- df1_common %>%
  arrange(padj) %>%
  slice(1:20) %>%
  pull(gene_id)

#making a combined matrix

combined_matrix <- cbind(
  A_vs_C = df1_common$log2FoldChange[match(top20_common, df1_common$gene_id)],
  A_vs_D = df2_common$log2FoldChange[match(top20_common, df2_common$gene_id)]
)

rownames(combined_matrix) <- top20_common

pheatmap(
  combined_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main= "combined heatmap",
  color = colorRampPalette(c("blue", "yellow", "pink"))(100)
)

#Significant gene lists or tables

# for A vs C
#upregulated genes
upreg_df1 <- df1 %>%
  filter(padj < 0.05, log2FoldChange > 0) %>%
  select(gene_id, log2FoldChange, pvalue, padj) %>%
  arrange(padj)

#downregulated genes
downreg_df1 <- df1 %>%
  filter(padj < 0.05, log2FoldChange < 0) %>%
  select(gene_id, log2FoldChange, pvalue, padj) %>%
  arrange(padj)

head(upreg_df1)
print(upreg_df1)
head(downreg_df1)
print(downreg_df1)

# for A vs D
#upregulated genes
upreg_df2 <- df2 %>%
  filter(p,adj < 0.05, log2FoldChange > 0) %>%
  select(gene_id, log2FoldChange, pvalue, padj) %>%
  arrange(padj)

#downregulated genes
downreg_df2 <- df2 %>%
  filter(padj < 0.05, log2FoldChange < 0) %>%
  select(gene_id, log2FoldChange, pvalue, padj) %>%
  arrange(padj)

head(upreg_df2)
print(upreg_df2)
head(downreg_df2)
print(downreg_df2)

#aditional analysis
#clustering of gene expression patterns
#a hierarchal clustering heatmap was already generated in the heatmap section that had the top 20 differentially expressed genes common between both datasets
