library(dplyr)
library(ggplot2)
plot_parent_weight <- function(){
  samples <- make('samples')
  
  mw <- samples %>%
    filter(!is.na(molwt_highlow)) %>%
    group_by(unique_id, molwt_highlow) %>%
    summarize(totals = sum(RESULT), counts = n()) %>%
    rename(variable = molwt_highlow)
  
  parent <- samples %>%
    filter(!is.na(parentAlkyl)) %>%
    group_by(unique_id, parentAlkyl) %>%
    summarize(totals = sum(RESULT), counts = n()) %>%
    rename(variable = parentAlkyl)
  
  all <- bind_rows(mw, parent)
  
  all$variable <- factor(all$variable, levels = c('parent', 'alkyl', 'LMW', 'HMW'))
  all.text <- ungroup(all) %>%
    select(variable, counts) %>%
    distinct() %>%
    mutate(totals = 1000000)
  
  p <- ggplot(all, aes(x = variable, y = totals)) +
    scale_y_continuous(trans = 'log10') +
    geom_boxplot() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = "Compound Type", y = "Sample Concentration (ppb)") +
    geom_text(data = all.text, aes(x = variable, y = totals, label = paste0("n = ", counts)))
  
  ggsave('9_parent_weight/doc/parent_mw_pah.png', p)
  
}
