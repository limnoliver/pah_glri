conc_by_site <- function(sample_dat, fig_file_path) {
  df <- sample_dat
  #df$unique_id <- paste0(df$State, "-", df$STAT_ID)
  
  df <- filter(df, PARAM_SYNONYM %in% "Priority Pollutant PAH") %>%
    select(unique_id, RESULT, State, Watershed, STAID, Lake, PARAM_SYNONYM) %>%
    group_by(unique_id, State, Watershed, STAID, Lake, PARAM_SYNONYM) %>%
    summarize(RESULT = mean(RESULT))
  
  df$STAID <- as.numeric(df$STAID)
  
 
  dfother <- filter(df, Lake != 'Michigan')
  dfmich <- filter(df, Lake == 'Michigan')
  # p <- pah::plot_pah(pah_dat = df, conc_column = 'RESULT', sample_id_column = 'unique_id',
  #                    compound_column = 'PARAM_SYNONYM', compound_plot = "Priority Pollutant PAH",
  #                    color_column = 'Lake', group_column = 'Watershed', order_column = 'STAID', conc_units = 'ppb')
  # 
  # p <- pah::plot_pah(pah_dat = df, conc_column = 'RESULT', sample_id_column = 'unique_id',
  #                    compound_column = 'PARAM_SYNONYM', compound_plot = "Priority Pollutant PAH",
  #                    color_column = 'Lake', order_column = 'STAID', conc_units = 'ppb')
  # 
  # p_other <- pah::plot_pah(pah_dat = dfother, conc_column = 'RESULT', sample_id_column = 'unique_id',
  #                          compound_column = 'PARAM_SYNONYM', compound_plot = "Priority Pollutant PAH",
  #                          color_column = 'Lake', group_column = 'Watershed', order_column = 'STAID', conc_units = 'ppb')
  # 
  # p_mich <- pah::plot_pah(pah_dat = dfmich, conc_column = 'RESULT', sample_id_column = 'unique_id',
  #                         compound_column = 'PARAM_SYNONYM', compound_plot = "Priority Pollutant PAH",
  #                         color_column = 'Lake', group_column = 'Watershed', order_column = 'STAID', 
  #                         conc_units = 'ppb')
  # 
  # barchart_theme <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  #                         legend.position = "top", panel.spacing=unit(5,"pt"), strip.background = element_blank(),
  #                         strip.text = element_text(margin = margin(4,10,4,10, unit="pt")))
  p_alt_other <- plot_pah_alt(pah_dat = dfother, all_dat = df, conc_column = 'RESULT', sample_id_column = 'unique_id',
                              compound_column = 'PARAM_SYNONYM', compound_plot = "Priority Pollutant PAH",
                              color_column = 'Lake', group_column = 'Watershed', group2_column = 'Lake', 
                              order_column = 'STAID', conc_units = 'ppb', guides = FALSE) +
    theme(plot.margin = margin(0.1,0,0,0, unit = 'cm'))
  
  p_alt_mich <- plot_pah_alt(pah_dat = dfmich, all_dat = df, conc_column = 'RESULT', sample_id_column = 'unique_id',
                              compound_column = 'PARAM_SYNONYM', compound_plot = "Priority Pollutant PAH",
                              color_column = 'Lake', group_column = 'Watershed', guides = FALSE,
                             group2_column = 'Lake', order_column = 'STAID', conc_units = 'ppb') +
    theme(plot.margin = margin(0,0,0,0, unit = 'cm'))
    
  
  
  dfall <- df
  dfall$Lake <- factor(dfall$Lake, levels = c('Erie', 'Huron', 'Ontario', 'Superior', 'Michigan'))
  
  p_alt_all <- plot_pah_alt(pah_dat = dfall, all_dat = df, conc_column = 'RESULT', sample_id_column = 'unique_id',
                            compound_column = 'PARAM_SYNONYM', compound_plot = "Priority Pollutant PAH",
                            color_column = 'Lake', group_column = 'Watershed', group2_column = 'Lake', 
                            order_column = 'STAID', conc_units = 'ppb') +
    theme(legend.position = "right", legend.box.margin = margin(0.2, 0.1, 0.1, 0.3, unit = 'cm'),
          legend.key.size = unit(0.4, "cm"))
  
  #library(cowplot)
  label1 <- "PEC = 22800 ppb"
  label2 <- "TEC = 1610 ppb"
  pother <- cowplot::ggdraw(p_alt_other) +
    cowplot::draw_label(label1, x = 0.09, y = 0.93, size = 10, hjust = 0) +
    cowplot::draw_label(label2, x = 0.09, y = 0.84, size = 10, hjust = 0)
  
  pmich <- cowplot::ggdraw(p_alt_mich) +
    cowplot::draw_label(label1, x = 0.105, y = 0.965, size = 10, hjust = 0) +
    cowplot::draw_label(label2, x = 0.105, y = 0.895, size = 10, hjust = 0)
  
  plegend <- cowplot::get_legend(p_alt_all)
  
  p1 <- cowplot::plot_grid(pmich, plegend, rel_widths = c(1, 0.14), align = 'h', axis = 't')
  p2 <- cowplot::plot_grid(pother, p1, align = 'v', ncol = 1, axis = 'l',
                     rel_heights = c(0.87, 1))
  
  ggplot2::ggsave(fig_file_path, 
         p2, height = 6, width = 7.25)
  
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
  
  toc_dat <- sample_dat %>%
    filter(PARAM_SYNONYM == "TOC") %>%
    mutate(f_TOC = RESULT/100) %>%
    select(unique_id, f_TOC)
  
  esbtu_dat <- sample_dat %>%
    filter(!is.na(coc_pah_fcv)) %>%
    left_join(toc_dat, by = 'unique_id') %>%
    mutate(conc_ug_g = round(RESULT/f_TOC, 0)/1000) %>%
    mutate(esbtu = conc_ug_g/coc_pah_fcv) %>%
    group_by(unique_id) %>%
    summarize(n_esbtu = n(), 
              sum_esbtu = sum(esbtu, na.rm = T))
  
  site_results <- left_join(pec_tec, esbtu_dat, by = 'unique_id') %>%
    rename(epa_priority16 = RESULT)
  
  site_results_summary <- data.frame(unique_id = c("TEC", 'PEC', 'ESBTU'), 
                                mean_EPApriority16_conc = rep(mean(site_results$epa_priority16), 3),
                                n_sites = rep(nrow(site_results), 3),
                                mean_ratio = c(mean(site_results$tec_ratio), mean(site_results$pec_ratio), mean(site_results$sum_esbtu)),
                                median_ratio = c(median(site_results$tec_ratio), median(site_results$pec_ratio), median(site_results$sum_esbtu)),
                                sd_ratio = c(sd(site_results$tec_ratio), sd(site_results$pec_ratio), sd(site_results$sum_esbtu)), 
                                min_ratio = c(min(site_results$tec_ratio), min(site_results$pec_ratio), min(site_results$sum_esbtu)),
                                max_ratio = c(max(site_results$tec_ratio), max(site_results$pec_ratio), max(site_results$sum_esbtu)))
  
  out <- list(site_results, site_results_summary)
  names(out) <- c('results_bysample', 'results_summary')
  return(out)
}