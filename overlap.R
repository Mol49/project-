

# import filtered datsets for each cell line (bigfc and sig)
#
# Jurkat 
Jurkat <- read.csv("results/bigfc_sig0.05_Jurkat.csv")
#
# KMBC2
KMBC2 <- read.csv("results/bigfc_sig0.05_KMBC2.csv")
#
# M231
M231 <- read.csv("results/bigfc_sig0.05_M231.csv")
#
# MCF&
MCF7 <- read.csv("results/bigfc_sig0.05_MCF7.csv")
#
###############################################################################
#
# Extract all gene symbols for each cell line
#
# Jurkat
Jurkat_genes <- Jurkat$symbol
#
# KMBC2
KMBC2_genes <- KMBC2$symbol
#
# M231
M231_genes <- M231$symbol
#
# MCF7
MCF7_genes <- MCF7$symbol
#
################################################################################
#
# Find overlapping genes in all cell lines
#
overlap_all <- Reduce(intersect, list(Jurkat_genes, KMBC2_genes, M231_genes, MCF7_genes))
overlap_all
# 0 genes overlap in all 4 cell lines
#
# find overlaps of 2 groups
#
pair_Jurkat_KMBC2 <- intersect(Jurkat_genes, KMBC2_genes)
pair_Jurkat_KMBC2 
# "TSC22D3"
pair_Jurkat_M231  <- intersect(Jurkat_genes, M231_genes)
pair_Jurkat_M231
# "NEBL"    "TSC22D3"
pair_Jurkat_MCF7  <- intersect(Jurkat_genes, MCF7_genes)
pair_Jurkat_MCF7
# 0
#
pair_KMBC2_M231   <- intersect(KMBC2_genes, M231_genes)
pair_KMBC2_M231
# "FKBP5"   "RIPOR2"  "TFCP2L1" "EDN2"    "VSTM2L"  "ADAMTS8" "TSC22D3" "TMEM63C"
# "SCNN1G" "SAA1"    "PER1"    "ALOX15B" "ASS1P1"  NA        "TRNP1"  
pair_KMBC2_MCF7   <- intersect(KMBC2_genes, MCF7_genes)
pair_KMBC2_MCF7
# "FKBP5" "SGK1"  "PER1"  NA 
#
pair_M231_MCF7    <- intersect(M231_genes, MCF7_genes)
pair_M231_MCF7
# "PDK4"   "FKBP5"  "ZBTB16" "PER1" 
#
# Find overlaps between 3 groups
# 
overlap_JKM <- Reduce(intersect, list(Jurkat_genes, KMBC2_genes, M231_genes))
overlap_JKM
# "TSC22D3"
overlap_JKM7 <- Reduce(intersect, list(Jurkat_genes, KMBC2_genes, MCF7_genes))
overlap_JKM7
# 0
overlap_JM7M <- Reduce(intersect, list(Jurkat_genes, MCF7_genes, M231_genes))
overlap_JM7M
# 0
overlap_KM7M <- Reduce(intersect, list(KMBC2_genes, MCF7_genes, M231_genes))
overlap_KM7M
# "FKBP5" "PER1" 
# 
###############################################################################
#
# make a 4 way venn plot
#
library(VennDiagram)
# 
# Remove NA values from gene lists
Jurkat_genes <- na.omit(Jurkat_genes)
KMBC2_genes  <- na.omit(KMBC2_genes)
M231_genes   <- na.omit(M231_genes)
MCF7_genes   <- na.omit(MCF7_genes)
# 
venn.plot <- venn.diagram(
  x = list(
    Jurkat = Jurkat_genes,
    KMBC2  = KMBC2_genes,
    M231   = M231_genes,
    MCF7   = MCF7_genes
  ),
  filename = "4way_venn_diagram.png",
  imagetype = "png",
  height = 2400,
  width = 2400,
  resolution = 300,
  col = "black",
  fill = c("red", "blue", "green", "yellow"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("red", "blue", "green", "yellow")
)
# 
###############################################################################
# 
# not many significant values in Jurkat cell line so try with KMBC2, MCF7 and 
# M231 :
# 
# Find all overlapping genes
overlapKMBC2.M231 <- KMBC2_genes[KMBC2_genes %in% M231_genes]
overlap.all <- MCF7_genes[MCF7_genes %in% overlapKMBC2.M231]
#
overlap.all
# there are 2 genes that overlap in all 3 datasets, FKBP5 and PER1
#
# overlap between KMBC2 and M231
#
overlapKMBC2.M231
# there are 15 genes that overlap, "FKBP5"   "RIPOR2"  "TFCP2L1" "EDN2"  "VSTM2L" 
# "ADAMTS8" "TSC22D3" "TMEM63C" "SCNN1G" "SAA1"    "PER1"    "ALOX15B" "ASS1P1" 
# "TRNP1"   "SAA1"   
#
# overlap between KMBC2 and MCF7
#
overlapKMBC2.MCF7 <- KMBC2_genes[KMBC2_genes %in% MCF7_genes]
overlapKMBC2.MCF7
# there are 3 genes that overlap, "FKBP5" "SGK1"  "PER1" 
#
# overlap between M231 and MCF7
#
overlapM231.MCF7 <- M231_genes[M231_genes %in% MCF7_genes]
overlapM231.MCF7
# there are 4 genes that overlap, "PDK4"   "FKBP5"  "ZBTB16" "PER1"  
#
################################################################################
# 
# make a 3 way venn plot
#
venn.plot.3way <- venn.diagram(
  x = list(
    KMBC2 = KMBC2_genes,
    M231  = M231_genes,
    MCF7  = MCF7_genes
  ),
  filename = "3way_venn_diagram.png",
  imagetype = "png",
  height = 2400,
  width = 2400,
  resolution = 300,
  col = "black",
  fill = c("blue", "green", "yellow"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("blue", "green", "yellow")
)
