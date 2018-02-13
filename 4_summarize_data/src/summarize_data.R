conc_by_site <- function(sample_dat, fig_file_path) {
  df <- sample_dat
  df$unique_id <- paste0(df$State, "-", df$STAT_ID)
  
  df <- filter(df, PARAM_SYNONYM %in% "Priority Pollutant PAH") %>%
    select(unique_id, RESULT, State, Watershed, Lake, PARAM_SYNONYM) %>%
    group_by(unique_id, State, Watershed, Lake, PARAM_SYNONYM) %>%
    summarize(RESULT = mean(RESULT))
  
  p <- pah::plot_pah(pah_dat = df, conc_column = 'RESULT', sample_id_column = 'unique_id',
                compound_column = 'PARAM_SYNONYM', compound_plot = "Priority Pollutant PAH",
                color_column = 'Lake', group_column = 'Watershed', conc_units = 'ppb')
  
  ggplot2::ggsave(fig_file_path, 
         p, height = 6, width = 14)
  
  # had to add "PEC" and "TEC" labels manually because
  # facets weren't wide enough to accomodate the text.
}