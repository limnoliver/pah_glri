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

