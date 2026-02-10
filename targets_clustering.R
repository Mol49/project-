# cluster analysis 2
#--------------------
#
hclust <- readRDS("results/hclust.RDS")
#
# cut hclust at height = 40, about halfway up the dendrogram.
clusters_h40 <- cutree(hclust, h = 40)
#
# clusters is a named vector with protiens and values as cluster IDs
clusters_h40[1:10]
# 
#           A0A0B4J2F2;P57059           A0A1B0GTU1;O75152           A0A1B0GUS4;P68036 
#                           1                           2                           3 
# A0A5B9;P0DSE2;P0DTU4;P01850                      A0AV02                      A0AV96 
#                           4                           1                           5 
#                      A0AVF1                      A0AVK6                      A0AVT1 
#                           3                           6                           7 
#                      A0FGR8 
#                           3 
# 
## See how many clusters height = 40 produces 
range(clusters_h40)
# = 21
# 
# Glucocorticoid receptor 
#---------------------------
# 
# the gene symbol for GR is NR3C1. Whats clustered with this?
#  NR3C1 is P04150
#
clusters_k1000 <- cutree(hclust, k = 1000)
clusters_k1000["P04150"]
# 
# GR (NR3C1) is in cluster 104
#
length(clusters_k1000[clusters_k1000 == 104])
# there are 15 proteins in this cluster, 
#
# what else is in this cluster?
c_k1000_104 <- names(clusters_k1000[clusters_k1000 == 104])
c_k1000_104
#
# "O00410"               "O14647"              "P04150"               "P34896"              
# "P36578"               "P52292"              "Q13951"               "Q14152"              
# "Q15542"               "Q8N0Y7;P15259;P18669""Q8NBZ7"               "Q92673"              
# "Q96AG4"               "Q9HAU0"              "Q9UNS2" 
#
# some seem relevant, such as importin-5 (O00410), which is involved in nuclear import
# or the TFIID subunit 5 (Q15542) which is involved in transcription initiation
# but others don't seem at all relevant.
#
clusters_h40 <- cutree(hclust, h = 40)
clusters_h40["P04150"] # cluster 17
length(clusters_h40[clusters_h40 == clusters_h40["P04150"]]) # 405 proteins in cluster 17
#
# ADAMTS8
#---------
#  ADAMTS8 is Q9UP79
#
clusters_k1000 <- cutree(hclust, k = 1000)
clusters_k1000["Q9UP79"]
# 
# ADAMST8 is in cluster NA

# RIPOR2
#--------
# 
# RIPOR2 is Q9Y4F9
# 
clusters_k1000 <- cutree(hclust, k = 1000)
clusters_k1000["Q9Y4F9"]
#
#  cluster NA
#
# PER1
#-----
#
# PER1 is O15534

clusters_k1000 <- cutree(hclust, k = 1000)
clusters_k1000["O15534"]
#
# PER1 is in cluster 5
#
length(clusters_k1000[clusters_k1000 == 5])
# there are 650 proteins in cluster 5
#





clusters_h40 <- cutree(hclust, h = 40)
clusters_h40["O15534"]
length(clusters_h40[clusters_h40 == clusters_h40["O15534"]])

# [1] 1113


#####################################
#
# overlap KMBC2 + M231
  #"FKBP5"   "RIPOR2"  "TFCP2L1" "EDN2"  "VSTM2L" "ADAMTS8" "TSC22D3"
# "TMEM63C" "SCNN1G" "SAA1"    "PER1"    "ALOX15B" "ASS1P1" "TRNP1"   "SAA1" 

# assign components to object overlap
# GO annotate says the following :
#"FKBP5"- Q13451
#"RIPOR2" - Q9Y4F9
#"TFCP2L1" -Q9NZI6
#"EDN2"  -P20800
#"VSTM2L" -Q96N03
#"ADAMTS8"-Q9UP79
#"TSC22D3"-Q99576
#"TMEM63C" -Q9P1W3
#"SCNN1G" -P51170
#"SAA1"  -P0DJI8
#"PER1" -O15534
#"ALOX15B" -O15296
#"ASS1P1" -
#"TRNP1"  -Q6NT89

overlap<- c("Q13451","Q9Y4F9","Q9NZI6","P20800","Q96N03","Q9UP79","Q99576","Q9P1W3","P51170","P0DJI8","O15534","O15296","Q6NT89")
clusters_h40[overlap]

# Q13451   <NA> Q9NZI6   <NA>   <NA>   <NA>   <NA> Q9P1W3   <NA>   <NA> O15534 
#      5     NA      1     NA     NA     NA     NA      1     NA     NA      1 
# O15296   <NA> 
#      1     NA 
#
# 4 of these are in cluster 1, ("TFCP2L1" 'TMEM63C' 'PER1' "ALOX15B")
# how many proteins are in cluster 1?
length(clusters_h40[clusters_h40 == 1])
# 1113
#
# if we cluster at h = 20, are these still in the same cluster?
clusters_h20 <- cutree(hclust, h = 20)
range(clusters_h20)
# there are now 96 clusters
#
clusters_h20[overlap]
length(clusters_h20[clusters_h20 == 1])
#  still in the same cluster and there are 779 proteins in this cluster

# h15
clusters_h15 <- cutree(hclust, h = 15)
range(clusters_h15)
# there are now 177 clusters
# 
clusters_h15[overlap]
#Q13451   <NA> Q9NZI6   <NA>   <NA>   <NA>   <NA> Q9P1W3   <NA>   <NA> O15534 
#   141     NA      5     NA     NA     NA     NA      5     NA     NA      5 
#O15296   <NA> 
#     5     NA 

length(clusters_h15[clusters_h15 == 5])









# what about slightly fewer clusters
clusters_h15 <- cutree(hclust, h = 15)
range(clusters_h15)
# 177 clusters
clusters_h15[AP1]
# still clustered!
length(clusters_h15[clusters_h15 == 120])
# 37 proteins in this cluster

# h = 11
clusters_h11 <- cutree(hclust, h = 11)
range(clusters_h11)
# 323 clusters
clusters_h11[AP1]
# still clustered
length(clusters_h11[clusters_h11 == 163])
# still 37

# so JUN and FOSL2 are clustered until h = 11, and their smallest joint cluster is a cluster of 37 proteins


