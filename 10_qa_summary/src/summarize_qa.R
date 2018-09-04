assess_dups <- function(qa_df = duplicates) {
  
  qa <- qa_df %>%
    mutate(unique_id = paste0(State, "-", STAT_ID)) %>%
    filter(!is.na(Parameter)) %>%
    select(unique_id, FIELD_QC_CODE, PARAM_SYNONYM, RESULT, molwt)
  
  dup <- filter(qa, FIELD_QC_CODE != "FB") %>%
    spread(key = FIELD_QC_CODE, value = RESULT) %>%
    mutate(perc_diff = round(abs((SA-DU)/ifelse(SA>=DU, SA, DU))*100, 2),
           one_bdl = (SA==0 & DU > 0)|(SA>0 & DU == 0))
  
  count_bdl_off <- dup %>%
    group_by(PARAM_SYNONYM) %>%
    summarize(n_both_bdl = length(which(is.nan(perc_diff))),
              n_one_bdl = length(which(one_bdl)), 
              n_duplicates = n())
  
  dup.stats <- filter(dup, one_bdl == FALSE) %>%
    group_by(PARAM_SYNONYM) %>%
    summarize(mean = mean(perc_diff, na.rm = T),
              median = median(perc_diff, na.rm = T),
              stdev = sd(perc_diff, na.rm = T), 
              min = min(perc_diff, na.rm = T), 
              max = max(perc_diff, na.rm = T))
  
  all.pah.means <- dup.stats$mean
  dup.stats[nrow(dup.stats) +1, 1] <- c("All PAH compounds")
  dup.stats[nrow(dup.stats), 2:6] <- c(mean(all.pah.means), 
                                       median(all.pah.means),
                                       sd(all.pah.means),
                                       min(all.pah.means),
                                       max(all.pah.means))
  
  dup.stats <- left_join(dup.stats, count_bdl_off)
  
  dup.stats[nrow(dup.stats), 7:9] <- c(sum(dup.stats$n_both_bdl, na.rm = T), sum(dup.stats$n_one_bdl, na.rm = T), sum(dup.stats$n_duplicates, na.rm = T))
  
  return(dup.stats)
}

assess_blanks <- function(qa_df) {
  qa <- qa_df %>%
    mutate(unique_id = paste0(State, "-", STAT_ID)) %>%
    filter(!is.na(Parameter)) %>%
    select(unique_id, FIELD_QC_CODE, PARAM_SYNONYM, RESULT, molwt)
  
  blank <- filter(qa, FIELD_QC_CODE != "DU") %>%
    spread(key = FIELD_QC_CODE, value = RESULT) %>%
    mutate(blank_perc_sample = round((FB/SA)*100, 2),
           FB_bdl = (FB == 0 & SA > 0),
           both_bdl = (FB == 0 & SA == 0),
           SA_bdl = (FB > 0 & SA == 0),
           both_adl = (FB > 0 & SA > 0))
  
  blank_bdl_counts <- group_by(blank, PARAM_SYNONYM) %>%
    summarize(n_field_blanks = n(), 
              n_field_blanks_bdl = length(which(FB_bdl == TRUE)),
              n_both_bdl = length(which(both_bdl == TRUE)),
              n_samples_bdl = length(which(SA_bdl == TRUE)),
              n_both_adl = length(which(both_adl == TRUE)))
  
  blank_perc_sample <- filter(blank, both_adl == TRUE) %>%
    group_by(PARAM_SYNONYM) %>%
    summarize(mean = mean(blank_perc_sample),
              median = median(blank_perc_sample),
              stdev = sd(blank_perc_sample),
              min = min(blank_perc_sample), 
              max = max(blank_perc_sample))
  
  blank_qa <- left_join(blank_bdl_counts, blank_perc_sample)
  return(blank_qa)
}

battelle_5433_qa_plots <- function() {
  
  ## lab control samples or reagent spikes
  qa_df <- make('pct_rec_labcontrol')
  
  qa_df <- filter(qa_df, !is.na(pct_recovery)) %>%
    filter(compound != "Biphenyl")
  
  wt_order <- select(samples, PARAM_SYNONYM, molwt) %>%
    distinct() %>%
    filter(PARAM_SYNONYM %in% unique(qa_df$compound)) %>%
    filter(PARAM_SYNONYM != "Biphenyl") %>%
    arrange(molwt)
  
  qa_df$compound <- factor(qa_df$compound, levels = wt_order$PARAM_SYNONYM)
  
  p <- ggplot(qa_df, aes(x = compound, y = pct_recovery)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "PAH compounds (low to high mol wt)", y = "% recovery", title = "Battelle lab control samples (n = 8)")
  ggsave('10_qa_summary/doc/battelle_lcs_pctrecovery.png', p, height = 4, width = 6)
  #################
  pct_rec_mspikes_5433 <- make('pct_rec_mspikes_5433')
  pct_rec_mspikes_5433$Parameter <- gsub('\\[', '\\(', x = pct_rec_mspikes_5433$Parameter)
  pct_rec_mspikes_5433$Parameter <- gsub('\\]', '\\)', x = pct_rec_mspikes_5433$Parameter)
  
  params_keep <- pct_rec_mspikes_5433$Parameter[which(pct_rec_mspikes_5433$Parameter %in% unique(qa_df$compound))]
  
  s_5433 <- filter(pct_rec_mspikes_5433, Parameter %in% params_keep) %>%
    mutate(source = 'NWQL_5433')
  
  #names(s_5433)[1:5] <- c('compound', 'min', 'mean', 'max', 'n')
  wt_order <- select(samples, PARAM_SYNONYM, molwt) %>%
    distinct() %>%
    filter(PARAM_SYNONYM %in% unique(s_5433$Parameter)) %>%
    filter(PARAM_SYNONYM != "Biphenyl") %>%
    arrange(molwt)
  
  s_5433$Parameter <- factor(s_5433$Parameter, levels = wt_order$PARAM_SYNONYM)
  
  p <- ggplot(s_5433, aes(x = Parameter, y = Median)) +
    geom_boxplot(aes(lower = First.Quartile, upper = Third.Quartile, middle = Median, ymin = Minimum, ymax = Maximum), stat = 'identity') +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "PAH compounds (low to high mol wt)", y = "% recovery", title = "NWQL schedule 5433 reagent spikes (n = 38)")
  ggsave('10_qa_summary/doc/nwql5433_reagspike_pctrecovery.png', p, height = 4, width = 6)
  
  ################ surrogates ####
  sur_5433 <- make('pct_rec_surrogates_5433')
  
  sur_5433$parameter <- gsub(",\\ssurrogate.*", "", sur_5433$parameter_nm)
  
  p <- ggplot(sur_5433, aes(x = parameter, y = result_va)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = "% recovery", title = "NWQL (5433) surrogate recovery (n = 70)")
  ggsave('10_qa_summary/doc/NWQL5433_surrogates_pctrecovery.png', p, height = 4, width = 4)
  
  
  sur_battelle <- make('pct_rec_surrogates')
  
   p <- ggplot(sur_battelle, aes(x = PARAM_SYNONYM, y = RESULT)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = '', y = "% recovery", title = "Battelle surrogate recovery (n = 70)")
  
  ggsave('10_qa_summary/doc/battelle_surrogates_pctrecovery.png', p, height = 4, width = 4)
   
  ######### matrix spikes ######
  
  battelle_ms <- make('pct_rec_mspikes')
  
  battelle_ms <- filter(battelle_ms, !is.na(pct_recovery)) %>%
    filter(compound != "Biphenyl")
  
  wt_order <- select(samples, PARAM_SYNONYM, molwt) %>%
    distinct() %>%
    filter(PARAM_SYNONYM %in% unique(battelle_ms$compound)) %>%
    filter(PARAM_SYNONYM != "Biphenyl") %>%
    arrange(molwt)
  
  battelle_ms$compound <- factor(battelle_ms$compound, levels = wt_order$PARAM_SYNONYM)
  
  # add n values to boxplot
  give.n <- function(x){
    return(c(y = median(x), label = length(x)))
  } 
  
  p <- ggplot(battelle_ms, aes(x = compound, y = pct_recovery)) +
    geom_boxplot() +
    stat_summary(fun.data = give.n, geom = "text", vjust = -0.3, color = 'red', size = 3) +
    coord_cartesian(ylim = c(40, 125)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = 'PAH compounds (low to high mol wt)', y = "% recovery", title = "Battelle matrix spike recovery (n values in red)")
  ggsave('10_qa_summary/doc/battelle_mspikes_pctrecovery.png', p, height = 4, width = 6)
  
}
