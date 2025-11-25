#
#
#
#
# Packages:
#----------
#
library(tidyverse)
#
library(readr)
#
library(GenomicRanges)
#
library(conflicted)
#
library(BiocManager)
#
library(readxl)
#
library(ggrepel)
#
library(pheatmap)
# 
###############################################################################
# 
# Set conflicts
#---------------
#
conflicts_prefer(dplyr::filter)
#
conflicts_prefer(GenomicRanges::setdiff)
#
###############################################################################
#
# Importing Data:
#----------------
#
Jurkat <- read_csv("results/results_Jurkat.csv")
#
M231 <- read_csv("results/results_M231.csv")
#
KMBC2 <- read_csv("results/results_KMBC2.csv")
#
MCF7 <- read_csv("results/results_MCF7.csv")
#
#################################################################################
#
# Recap of the results for Jurkat
#
glimpse(Jurkat)
# 
# 32,764 rows (genes) in 8 coloumns (baseMean, LFC, lfcSE, stat, pvalue, padj, symbol)
# 
# create dataframe of the significant genes and write them to file.
# This is subset from the results file but will make it a little easier to
# examine genes of interest.

# create a dataframe of genes significant at 0.05 level and save as csv
Jurkat_sig0.05 <- Jurkat |> 
  filter(padj <= 0.05)
#
# 7 genes are significant at 0.05 level
#
write_csv(Jurkat_sig0.05,
          "results/Jurkat_sig0.05.csv")
#
###############################################################################
#
# Recap of the results for KMBC2
#
glimpse(KMBC2)
# 
# 32,764 rows (genes) in 8 coloumns (baseMean, LFC, lfcSE, stat, pvalue, padj, symbol)
# 
# create dataframe of the significant genes and write them to file.
# This is subset from the results file but will make it a little easier to
# examine genes of interest.

# create a dataframe of genes significant at 0.05 level and save as csv
KMBC2_sig0.05 <- KMBC2 |> 
  filter(padj <= 0.05)
#
# 1141 genes are significant at 0.05 level
#
write_csv(KMBC2_sig0.05,
          "results/KMBC2_sig0.05.csv")

###############################################################################
#
# Recap of the results for M231
#
glimpse(M231)
# 
# 32,764 rows (genes) in 8 coloumns (baseMean, LFC, lfcSE, stat, pvalue, padj, symbol)
# 
# create dataframe of the significant genes and write them to file.
# This is subset from the results file but will make it a little easier to
# examine genes of interest.

# create a dataframe of genes significant at 0.05 level and save as csv
M231_sig0.05 <- M231 |> 
  filter(padj <= 0.05)
#
# 677 genes are significant at 0.05 level
#
write_csv(M231_sig0.05,
          "results/M231_sig0.05.csv")

###############################################################################
#
# Recap of the results for MCF7
#
glimpse(MCF7)
# 
# 32,764 rows (genes) in 8 coloumns (baseMean, LFC, lfcSE, stat, pvalue, padj, symbol)
# 
# create dataframe of the significant genes and write them to file.
# This is subset from the results file but will make it a little easier to
# examine genes of interest.

# create a dataframe of genes significant at 0.05 level and save as csv
MCF7_sig0.05 <- MCF7 |> 
  filter(padj <= 0.05)
#
# 182 genes are significant at 0.05 level
#
write_csv(MCF7_sig0.05,
          "results/MCF7_sig0.05.csv")

###############################################################################
#
# Volcano plotting
#------------------
#
# Volcano plots visualize the results following differential expression analysis
# we plot the negative log of the adjusted p value against the log 
# fold change.
# By plotting the negative log of the adjusted p-value the values are spread out,
# and the most significant are at the top of the axis
#
################################################################################
#
# Volcano plot for Jurkat
#-------------------------
# 
# Add a column to the results dataframe that contains the -log10(padj).
# Log10 is used so that the values are spread out and more easily visualized
#
Jurkat <- Jurkat |> 
  mutate(log10_padj = -log10(padj))
#
# specify in results if gene is significant (0.05) or fold change high
#
Jurkat <- Jurkat |> 
  mutate(sig = padj <= 0.05,
         bigfc = abs(log2FoldChange) >= 2)

# create data set containing only big fc and significant genes(0.05)
bigfc_sig0.05_Jurkat <- Jurkat %>%
  filter(sig == TRUE, bigfc == TRUE)
#
# save this to results folder
write_csv(bigfc_sig0.05_Jurkat, file = "results/bigfc_sig0.05_Jurkat.csv")
# 
# create volcano plot with labels for significant genes and large fold changes
#
Jurkat |> 
  ggplot(aes(x = log2FoldChange, 
             y = log10_padj, 
             colour = interaction(sig, bigfc))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed") +
  geom_vline(xintercept = 2, 
             linetype = "dashed") +
  geom_vline(xintercept = -2, 
             linetype = "dashed") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_manual(values = c("gray", 
                                 "lightgreen",
                                 "gray30",
                                 "green3")) +
  geom_text_repel(data = Jurkat |> 
                    filter(bigfc == TRUE, sig == TRUE),
                  aes(label = symbol),
                  size = 3,
                  max.overlaps = 50) +
  theme_classic() +
  theme(legend.position = "none")
#
# Save this plot
#
ggsave("figures/labelledJurkatvolcano.png",
       height = 4.5, 
       width = 4.5,
       units = "in",
       device = "png")
#
###############################################################################
#
# Volcano plot for KMBC2
#-------------------------
# 
# Add a column to the results dataframe that contains the -log10(padj).
# Log10 is used so that the values are spread out and more easily visualized
#
KMBC2 <- KMBC2 |> 
  mutate(log10_padj = -log10(padj))
#
# specify in results if gene is significant (0.05) or fold change high
#
KMBC2 <- KMBC2 |> 
  mutate(sig = padj <= 0.05,
         bigfc = abs(log2FoldChange) >= 2)

# create data set containing only big fc and significant genes(0.05)
bigfc_sig0.05_KMBC2 <- KMBC2 %>%
  filter(sig == TRUE, bigfc == TRUE)
#
# save this to results folder
write_csv(bigfc_sig0.05_KMBC2, file = "results/bigfc_sig0.05_KMBC2.csv")
# 
# create volcano plot with labels for significant genes and large fold changes
#
KMBC2 |> 
  ggplot(aes(x = log2FoldChange, 
             y = log10_padj, 
             colour = interaction(sig, bigfc))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed") +
  geom_vline(xintercept = 2, 
             linetype = "dashed") +
  geom_vline(xintercept = -2, 
             linetype = "dashed") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_manual(values = c("gray", 
                                 "lightgreen",
                                 "gray30",
                                 "green3")) +
  geom_text_repel(data = KMBC2 |> 
                    filter(bigfc == TRUE, sig == TRUE),
                  aes(label = symbol),
                  size = 3,
                  max.overlaps = 50) +
  theme_classic() +
  theme(legend.position = "none")
#
# Save this plot
#
ggsave("figures/labelledKMBC2volcano.png",
       height = 4.5, 
       width = 4.5,
       units = "in",
       device = "png")
