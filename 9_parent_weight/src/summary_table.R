summary_table <- function(thresholds, parent_mw, location) {
  thresholds <- make('threshold_dat')
  parent_mw <- make('parent_mw_dat')
  
  thresh_dat <- thresholds$results_bysample %>%
    select(-n_esbtu)
  
  type_dat <- parent_mw$summarized_dat %>%
    select(unique_id, parent_alkyl, HMW_LMW)
  
  summary_table <- left_join(thresh_dat, type_dat) %>%
    mutate(sum_EPA16 = round(sum_EPA16, 0),
           tec_ratio = round(tec_ratio, 2),
           pec_ratio = round(pec_ratio, 2),
           sum_esbtu = round(sum_esbtu, 2),
           perc_toc = round(perc_toc, 1), 
           parent_alkyl = round(parent_alkyl, 1), 
           HMW_LMW = round(HMW_LMW, 1))
  
  write.csv(x = summary_table, file = location, row.names = F)
}