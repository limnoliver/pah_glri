make_ratios <- function(pah_dat = samples){
  #pah_dat <- make('samples')
  pah_dat$sample_id <- paste0(pah_dat$State, "-", pah_dat$STAT_ID)
  pah_dat <- select(pah_dat, sample_id, casrn, PARAM_SYNONYM, RESULT) %>%
    group_by(sample_id, PARAM_SYNONYM, casrn) %>%
    summarize(RESULT = mean(RESULT)) %>%
    ungroup()
  
  ratios <- calc_ratios(pah_dat, sample_column = 'sample_id', conc_column = 'RESULT')
  return(ratios)
}

ratio_plot <- function(filename) {
  p <- plot_ratios(ratios)
  ggsave(filename, p, height = 5, width = 16)
}
