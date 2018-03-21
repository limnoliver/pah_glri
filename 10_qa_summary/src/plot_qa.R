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
              n_one_bdl = length(which(one_bdl)))
  
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
  
  dup.stats[nrow(dup.stats), 7:8] <- c(sum(dup.stats$n_both_bdl, na.rm = T), sum(dup.stats$n_one_bdl, na.rm = T))
  
  return(dup.stats)
}

# assess blanks

assess_blanks <- function(qa_df = blanks) {
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
    summarize(n_FB_bdl = length(which(FB_bdl == TRUE)),
              n_both_bdl = length(which(both_bdl == TRUE)),
              n_SA_bdl = length(which(SA_bdl == TRUE)),
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
