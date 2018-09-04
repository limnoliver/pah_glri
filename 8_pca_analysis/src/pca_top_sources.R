get_pca_top_sources <- function(pca_dat) {
  
  pca_dat <- pca_dat[[3]]
  top_sources <- group_by(pca_dat, sample) %>%
    summarize(top_source = source[which.min(euc_dist)])
  
  source.names <- pah::sources
  
  top_sources <- left_join(top_sources, source.names, by = c('top_source' = 'source_abbrev'))
  return(top_sources)
  
}

get_pca_top_sources_simple <- function(pca_dat) {
  
  pca_dat <- pca_dat[[3]]
  top_sources <- group_by(pca_dat, sample) %>%
    summarize(top_source = source[which.min(euc_dist)])
  
  return(top_sources)
  
}