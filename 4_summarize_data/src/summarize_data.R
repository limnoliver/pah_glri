conc_by_site <- function(sample_dat, fig_file_path) {
  df <- sample_dat
  df$unique_id <- paste0(df$State, "-", df$STAT_ID)
  
  df <- filter(df, PARAM_SYNONYM %in% "Priority Pollutant PAH") %>%
    select(unique_id, RESULT, State, Watershed, STAID, Lake, PARAM_SYNONYM) %>%
    group_by(unique_id, State, Watershed, STAID, Lake, PARAM_SYNONYM) %>%
    summarize(RESULT = mean(RESULT))
  
  df$STAID <- as.numeric(df$STAID)
  
  p <- pah::plot_pah(pah_dat = df, conc_column = 'RESULT', sample_id_column = 'unique_id',
                     compound_column = 'PARAM_SYNONYM', compound_plot = "Priority Pollutant PAH",
                     color_column = 'Lake', group_column = 'Watershed', order_column = 'STAID', conc_units = 'ppb')
  #library(cowplot)
  label1 <- "PEC = 22800 ppb"
  label2 <- "TEC = 1610 ppb"
  p2 <- cowplot::ggdraw(p) +
    cowplot::draw_label(label1, x = 0.055, y = 0.815, size = 14, hjust = 0) +
    cowplot::draw_label(label2, x = 0.055, y = 0.71, size = 14, hjust = 0)
  
  ggplot2::ggsave(fig_file_path, 
         p2, height = 6, width = 14)
  
  # had to add "PEC" and "TEC" labels manually because
  # facets weren't wide enough to accomodate the text.
}

tox_thresholes <- function(sample_dat) {
  sample_dat <- make('samples')
  
  tec <- ifelse(conc_unit == 'ppb', 1610, 1.610)
  pec <- ifelse(conc_unit == 'ppb', 22800, 22.8)
  
  pec_tec <- filter(sample_dat, PARAM_SYNONYM == "Priority Pollutant PAH") %>%
    mutate(tec_ratio = RESULT/tec,
           pec_ratio = RESULT/pec) %>%
    select(unique_id, RESULT, tec_ratio, pec_ratio)
  
  pec_tec_summary <- data.frame(unique_id = c("all sites - TEC", 'all sites - PEC'), 
                                mean_priority16_conc = rep(mean(pec_tec$RESULT), 2),
                                mean_ratio = c(mean(pec_tec$tec_ratio), mean(pec_tec$pec_ratio)),
                                median_ratio = c(median(pec_tec$tec_ratio), median(pec_tec$pec_ratio)),
                                sd_ratio = c(sd(pec_tec$tec_ratio), sd(pec_tec$pec_ratio)), 
                                min_ratio = c(min(pec_tec$tec_ratio), min(pec_tec$pec_ratio)),
                                max_ratio = c(max(pec_tec$tec_ratio), max(pec_tec$pec_ratio)))
  esbtu_dat <- sample_dat %>%
    filter(!is.na(coc_pah_max))
}