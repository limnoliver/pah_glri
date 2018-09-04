min_sources <- function(prof_dat) {

  top_sources <- prof_dat$sum_chi2 %>%
    group_by(unique_id) %>%
    summarize(top_source = source[which.min(sum_chi2)])
   
  top_sources 
    
}