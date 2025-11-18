###################################################################
#
####################
#
# Continuing analysis of RIME/RNAseq data
# ----------------------------------------
# 
# Cluster protiens based on how similar correlation profiles are across
# all trasncripts.
# hclus() will be used to create a dendrogram as heatmap() is incompatible
# with a matrix this large.
# https://www.datacamp.com/tutorial/hierarchical-clustering-R?dc_referrer=https%3A%2F%2Fwww.google.com%2F 
# 
# Import previously generated matrix
#
cor_stats <- readRDS(file = "results/mat2.rds")
# 
mat <- cor_stats
# 
# use hclust() to cluster the rows so the matrix can be transposed
#
mat <- t(mat)
dim(mat)
# 
# check for NA values and remove rows/cols with NAs
#
sum(is.na(mat))
mat <- mat[, colSums(is.na(mat)) == 0]
mat <- mat[rowSums(is.na(mat)) == 0, ]
sum(is.na(mat))
dim(mat)
#
# mat is a correlation matrix, hclust() requires a distance object
dist <- dist(mat)
hclust <- hclust(dist, method = "average")
saveRDS(hclust, "results/hclust.rds")
#
#####################################################################
