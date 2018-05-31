random_top_sources <- data.frame(unique_id = pca_top$sample, 
                                 top_source = sample(pah::sources$source_abbrev, size = length(pca_top$sample), replace = T))

write.csv(random_top_sources, '12_top_sources_by_site/raw/mix_mod_top_sources.csv', row.names = F)
