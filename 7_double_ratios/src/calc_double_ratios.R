make_ratios <- function(pah_dat = samples){
  #pah_dat <- make('samples')
  #pah_dat$sample_id <- paste0(pah_dat$State, "-", pah_dat$STAT_ID)
  pah_dat <- select(pah_dat, unique_id, casrn, PARAM_SYNONYM, RESULT) %>%
    group_by(unique_id, PARAM_SYNONYM, casrn) %>%
    summarize(RESULT = mean(RESULT)) %>%
    ungroup()
  
  ratios <- calc_ratios(pah_dat, sample_column = 'unique_id', conc_column = 'RESULT')
  return(ratios)
}

ratio_plot <- function(filename, ratio_dat) {
  p <- plot_ratios(ratio_dat)
  ggsave(filename, p, height = 5, width = 16)
}

order_samples <- function(totals) {
  sample_order <- totals %>%
    arrange(Priority16)

  sample_order <- sample_order$unique_id
  return(sample_order)
}

get_ratio_top_sources <- function(ratio_dist) {
  temp <- ratio_dist$sample
  return(temp)
}
