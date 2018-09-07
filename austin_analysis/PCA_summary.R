# Script to summarize PCA results 

library(remake)
library(dplyr)
library(pah)

delete("pca_creosote")
pca <- make("pca_creosote")

dist <- pca$pca_distance

distSummary <- group_by(dist, source) %>%
  summarize(n=n(),
            median = median(euc_dist, na.rm=TRUE),
            mean = mean(euc_dist, na.rm=TRUE))

write.csv(distSummary, "C:/Users/akbaldwi/Documents/GLRI/PAHs/pah/pah_glri/austin_analysis/PCA_euclidean_distances_summary.csv", row.names = F)


dat <- pca$pca_dat

varianceSummary <- pca$pca_summary

