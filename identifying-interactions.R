########################################################
# 
#########
#
# Introduction:
#----------------
# RIME data includes protiens that interact with the glucocorticoid
# receptor (GR) compared to an IgG control. The RIME data shows
# log fold changes for each protien, DE_results.csv was initially 
# processed using fragpipe. 
# RNA sequencing data of cells treated with and without Dex, an
# activator of GR, is analysed using DESeq2. RNA-seq data is found
# in Dex_RNA-seq.rds. 
# In this script i will analyse the log 2 fold changes of the two data
# sets (Log2FC of RNA vs. Log2FC of protiens interactions with GR) to
# determine correlations between protien and transcripts when GR is 
# activated.
#
########################################################################
#
# Set up:
#--------
#
options(repos = c(CRAN = "http://cran.rstudio.com"))
#
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#
BiocManager::install("DESeq2")
#
library(tidyverse)
#
library(DESeq2)
#
#########################################################################
# 
# Import RIME data:
#------------------
# 
rime_data <- read_csv("data-raw/DE_results.csv")
#
# Keep only the columns showing the fold change values of the
# same cell line, which also appear in the RNA-seq data set
# (Jurkat, MCF7, KMBC2, M231).
#
rime_lfcs <- rime_data[,(names(rime_data) %in% 
                           c("Protein ID",
                             "Gene Name",
                             "Jurkat_IgG_vs_Jurkat_GR_log2 fold change",
                             "KMBC2_IgG_vs_KMBC2_GR_log2 fold change",
                             "MCF7_IgG_vs_MCF7_GR_log2 fold change",
                             "X231_IgG_vs_X231_GR_log2 fold change"))]
#
# rename LogFCs to simpler names 
# 
names(rime_lfcs)[names(rime_lfcs) == "Jurkat_IgG_vs_Jurkat_GR_log2 fold change"] <- "Jurkat"
names(rime_lfcs)[names(rime_lfcs) == "KMBC2_IgG_vs_KMBC2_GR_log2 fold change"] <- "KMBC2"
names(rime_lfcs)[names(rime_lfcs) == "MCF7_IgG_vs_MCF7_GR_log2 fold change"] <- "MCF7"
names(rime_lfcs)[names(rime_lfcs) == "X231_IgG_vs_X231_GR_log2 fold change"] <- "231"
#
# # remove the gene ID (make new dataframe as gene names 
# could be useful later on)
#
rime_lfcs_p <- rime_lfcs[,(names(rime_lfcs) %in%
                             c("Protein ID",
                               "Jurkat",
                               "KMBC2",
                               "MCF7",
                               "231"))]
#
# set protein IDs as row names
rime_lfcs_p <- as.data.frame(rime_lfcs_p) # no longer a tibble
rownames(rime_lfcs_p) <- rime_lfcs_p$`Protein ID`
rime_lfcs_p$`Protein ID` <- NULL
#
# give this data frame a simple name for later
#
rime <- rime_lfcs_p
#
#########################################################################
#
# Import RNA-seq data:
#---------------------
#
dds <- readRDS("data-raw/Dex_RNA-seq.rds")
#
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("condition", "cellline"))
#
# add columns to dds combining condition and cellline 
# ( so we can compare logFCs of individual cell lines)
#
colData(dds)$condition_cellline <- factor(paste0(dds$condition, "_", dds$cellline))
#
# update dds design to include new combined factors
#
design(dds) <- ~ condition_cellline
#
# filter dataset to remove low expression counts and run DESeq2
#
dds <- dds[rowSums(counts(dds)) > 1,]
dds <- DESeq(dds)
#
# compare DMSO vs Dex treatments in Jurkat cells and 
# output the results as a dataframe 
#
res_Jurkat <- results(dds, contrast = c("condition_cellline", "DMSO_Jurkat", "Dex_Jurkat"))
results_Jurkat <- as.data.frame(res_Jurkat)
#
# add column wiht the cell line to the results 
#
results_Jurkat$cellline <-c("Jurkat")
#
# delete all columns except for celllines, padj,
# Log2foldchange
#
results_Jurkat <- results_Jurkat |>
  select(cellline, padj, log2FoldChange)
#
# rename log2foldchange to specify this is the RNA LFC
#
names(results_Jurkat)[names(results_Jurkat) == "log2FoldChange"] <- "rna_lfc"
#
# repeat for other cell lines that appear in both datasets
# (MCF7, KMBC2, M231) 
#
# MCF7
res_MCF7 <- results(dds, contrast = c("condition_cellline", "DMSO_MCF7", "Dex_MCF7"))
results_MCF7 <- as.data.frame(res_MCF7)
results_MCF7$cellline <- c("MCF7")
results_MCF7 <- results_MCF7 |>
  select(cellline, padj, log2FoldChange)
names(results_MCF7)[names(results_MCF7) == "log2FoldChange"] <- "rna_lfc"
#
# KMBC2
res_KMBC2 <- results(dds, contrast = c("condition_cellline", "DMSO_KMBC2", "Dex_KMBC2"))
results_KMBC2 <- as.data.frame(res_KMBC2)
results_KMBC2$cellline <- c("KMBC2")
results_KMBC2 <- results_KMBC2 |>
  select(cellline, padj, log2FoldChange)
names(results_KMBC2)[names(results_KMBC2) == "log2FoldChange"] <- "rna_lfc"
#
# 231
res_231 <- results(dds, contrast = c("condition_cellline", "DMSO_M231", "Dex_M231"))
results_231 <- as.data.frame(res_231)
results_231$cellline <- c("231")
results_231 <- results_231 |>
  select(cellline, padj, log2FoldChange)
names(results_231)[names(results_231) == "log2FoldChange"] <- "rna_lfc"
#
# merge all RNA-seq datasets into one 
#
rnaseq_lfcs <- bind_rows(results_Jurkat, results_KMBC2, results_231, results_MCF7)
#
# add column containing the gene names 
#
rnaseq_lfcs$gene <- rownames(rnaseq_lfcs)
rownames(rnaseq_lfcs) <- NULL
#
# reorder columns
#
rnaseq_lfcs <- rnaseq_lfcs[, c("gene", "cellline", "rna_lfc")]
#
# ensembl IDs, containing its version, are used for gene names 
# remove the version from the ensembl IDs and 
# # convert to a wide format. each column is a different cell type
#
# now put into wide format
#
rnaseq_wide <- rnaseq_lfcs |>
  pivot_wider(names_from = cellline, values_from = rna_lfc)
#
rnaseq_lfcs$gene <- rownames(rnaseq_lfcs)
rownames(rnaseq_lfcs) <- NULL
#
# reorder columns
rnaseq_lfcs <- rnaseq_lfcs[, c("gene", "cellline", "rna_lfc")]
#
# put gene names back into the column name
#
rnaseq_wide <- as.data.frame(rnaseq_wide) # no longer a tibble
rownames(rnaseq_wide) <- rnaseq_wide$gene
rnaseq_wide$gene <- NULL
#
# convert data frame to a matrix so its numeric and give a
# simpler name 
#
rnaseq <- as.matrix(rnaseq_wide)
# 
###########################################################################
# 
# Correlation test:
#------------------
# 
# create a matrix of compare all protiens and transcripts 
#
mat1 <- matrix(ncol=nrow(rime), nrow=nrow(rnaseq))
#
colnames(mat1) <- rownames(rime)
rownames(mat1) <- rownames(rnaseq)
#
mat2 <- mat1
#
# create a for loop which does a Pearson test for each 
# combination of protien and transcript.
# 
for (xsamples in rownames(rime)) {
  
  # convert rime row to numeric vector
  x <- rime[xsamples, ]
  
  cor_mat <- apply(rnaseq, 1, function(y) { 
    test <- cor.test(x, y)
    c(cor = test$estimate, p.value = test$p.value)
  } )
  
  cor_mat <- t(cor_mat)
  
  # putting names here will slow up a bit, but protects against reordering.
  mat1[names(cor_mat[,'cor.cor']), xsamples] <- cor_mat[,'cor.cor']
  mat2[names(cor_mat[,'p.value']), xsamples] <- cor_mat[,'p.value']
  
  
}

# Run in viking and save results 
#
saveRDS(mat1, "results/mat1.rds")
saveRDS(mat2, "results/mat2.rds")
#
############################################################################











