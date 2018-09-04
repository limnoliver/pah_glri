pot_sources <- make('ratio_distance')
pot_sources <- sources$sample
pot_sources <- c(sources$top_source1, sources$top_source2, sources$top_source3)

random_top_sources <- data.frame(unique_id = pca_top$sample, 
                                 top_source = sample(pot_sources, size = length(pca_top$sample), replace = T))

write.csv(random_top_sources, '12_top_sources_by_site/raw/mix_mod_top_sources.csv', row.names = F)
