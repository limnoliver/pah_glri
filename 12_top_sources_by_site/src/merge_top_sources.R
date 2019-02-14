merge_top_sources <- function(profiles_top, pca_top, pmf_dat, sample_order) {
  
  
  order <- data.frame(sample = sample_order, rank = 1:length(sample_order), stringsAsFactors = F)
  
  pca_top <- rename(pca_top, pca_top_sources = top_source, pca_distance = top_distance)
  profiles_top <- rename(profiles_top, profiles_top_sources = top_source, profiles_distance = top_distance)
  parent_weight <- parent_mw_dat[[2]] %>%
    select(unique_id, parent_alkyl, HMW_LMW)
  
  
  top_sources <- left_join(pca_top, profiles_top, by = c('sample' = 'unique_id')) %>%
    left_join(parent_weight, by = c('sample' = 'unique_id'))
  
  unique_top_sources <- unique(c(as.character(top_sources$pca_top_sources), top_sources$profiles_top_sources))
  
  
  
  
  
  #my.cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(8, 'RdBu')))
  my_cols = qualitative_hcl(length(unique_top_sources), "Dark 3")
  my.cols <- data.frame(my_cols = qualitative_hcl(length(unique_top_sources), "Dark 3"),
                        sources = unique_top_sources, stringsAsFactors = FALSE)
  
  top_sources <- top_sources %>%
    mutate(`HMW:LMW` = ifelse(HMW_LMW > 1, 'pyro', 'petro'),
           `Parent:Alkyl` = ifelse(parent_alkyl > 1, 'pyro', 'petro')) %>%
    left_join(rename(my.cols, pca_top_sources = sources, pca_colors = my_cols)) %>%
    left_join(rename(my.cols, profiles_top_sources = sources, profiles_colors = my_cols)) %>%
    left_join(select(pmf_top, sample = site, pmf_top_sources = dominantSource)) %>%
    left_join(order) %>%
    rename(`PCA top source` = pca_top_sources,
           `Profiles top source` = profiles_top_sources) %>%
    arrange(sample)
  
  
  return(top_sources)
  
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
