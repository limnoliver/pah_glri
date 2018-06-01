merge_top_sources <- function(profile_top = , ratio_top = , pca_top = , mix_top =) {
  pca_top$top_source <- as.character(pca_top$top_source)
  
  top_sources <- rename(ratio_top,
                        unique_id = sample, 
                        ratio1 = top_source1,
                        ratio2 = top_source2,
                        ratio3 = top_source3) %>%
    left_join(rename(profile_top, profile = top_source)) %>%
    left_join(rename(pca_top, unique_id = sample, pca = top_source)) %>%
    left_join(rename(mix_top, mix = top_source))
  
  return(top_sources)
  
  # output list of unique sources
  all.top <- unique(c(top_sources$ratio1, top_sources$ratio2, top_sources$ratio3, top_sources$profile, top_sources$pca))
  write.csv(all.top, '12_top_sources_by_site/doc/unique_top_sources.csv', row.names = F)
}

calc_mass_rules <- function(top = merged_top_sources, mass_rules = percent_by_weight_bysample) {
  top_long <- gather(top, key = method, value = top_source, -unique_id)
  mass_long <- gather(mass_rules, key = source_name, value = mass_fraction, -unique_id)
  
  # read in source crosswalk data
  cross <- read.csv('12_top_sources_by_site/raw/source_crosswalk.csv', stringsAsFactors = F, na.strings = '')
  
  top_long <- left_join(top_long, cross, by = c('top_source' = 'source_abbreviation'))
  
  top_with_mass <- left_join(top_long, mass_long)
  top_with_mass$mass_fraction <- as.numeric(gsub('>', '', top_with_mass$mass_fraction))
  
  impossible <- which(top_with_mass$mass_fraction == 100)
  unlikely <- which(top_with_mass$mass_fraction >= 5)
  
  p <- ggplot(top_with_mass, aes(x = method, y = unique_id)) + 
    geom_tile(aes(fill = top_source), color = "white")
  
  pb <- ggplot_build(p)
  p2 <- p +
    #geom_segment(data = pb$data[[1]][impossible, ],
    #             aes(x = xmin, xend = xmax, y = ymin, yend = ymax), size = 1.5) +
    geom_segment(data = pb$data[[1]][unlikely, ],
                 aes(x = xmin, xend = xmax, y = ymin, yend = ymax)) +
    theme_bw() + 
    labs(x = "Method", y = "")
  
  ggsave('12_top_sources_by_site/doc/top_source_by_site_summary.png', p2, height = 10, width = 7)
    
}
