# calculate difference between source and sample profiles
library(dplyr)
library(rlang)
library(tidyr)

samples$sample_id <- paste0(samples$State, "-", samples$STAT_ID)

# remove all samples where there is a zero - meaning BDL

samples <- filter(samples, sourceProfile12 == TRUE)
samples.remove <- unique(samples$sample_id[samples$RESULT == 0])
samples <- filter(samples, !(sample_id %in% samples.remove))
samples <- select(samples, sample_id, RESULT, casrn, Lake, State, Watershed)

pah_profiler <- function(pah_dat, compound_column = 'casrn', sample_column = 'sample_id', 
                         conc_column = 'RESULT', sources = source_profiles) {
  # make column names dplyr-ready
  quo_compound_column <- sym(compound_column)
  quo_conc_column <- sym(conc_column)
  quo_sample_column <- sym(sample_column)
  
  # pull out all 12 compounds
  profile_compounds <- select(sources, !!quo_compound_column)
  
  # filter user samples to include only those in the source profiles
  samp.prof <- filter(pah_dat, (!!quo_compound_column) %in% profile_compounds[[1]])
  
  # group by sample - calculate sum total and whether or not there are any samples = 0
  samp.prof.bysample <- group_by(samp.prof, !!quo_sample_column) %>%
    summarize(total_pah = sum(!!quo_conc_column))
  
  samp.prof <- left_join(samp.prof, samp.prof.bysample) %>%
    mutate(prop_conc = (!!quo_conc_column)/total_pah)
  
  # merge in source compound info
  all.profs <- full_join(samp.prof, sources, by = compound_column) %>%
    select(-RESULT) %>%
    gather(key = source, value = source_prop_conc, -(1:8), -Compound, -Abbreviation, -pcode, -molwt)
  
  # calculate the chi squared difference
  all.profs <- mutate(all.profs, chi2 = (abs(prop_conc - source_prop_conc)^2)/((prop_conc + source_prop_conc)/2))
  
  # sum over the sources by sample
  all.summary <- group_by(all.profs, sample_id, source) %>%
    summarize(sum_chi2 = sum(chi2))
  
  # plot source chi squareds across all samples
  # first show dropping creosote
  
  # order sources by median value
  creosote <- unique(grep('creosote', all.summary$source, ignore.case = T, value = T))
  all.no.creosote <- filter(all.summary, !(source %in% creosote))
  order.vals <- all.no.creosote %>%
    group_by(source) %>%
    summarize(med = median(sum_chi2)) %>%
    arrange(med)
  
  all.no.creosote$source <- factor(all.no.creosote$source, levels = order.vals$source)
  
  p <- ggplot(all.no.creosote, aes(x = source, y = sum_chi2)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = '', y = 'Sum Chi2')
  
  ggsave('6_profile_data/doc/chi2_allsites_nocreosote.png', p)
  
  
  compound.drop <- unique(all.profs$casrn[is.na(all.profs$source_prop_conc)])
  creosote.summary <- filter(all.profs, casrn != compound.drop) %>%
    group_by(sample_id, source) %>%
    summarize(sum_chi2 = sum(chi2))
  
  order.vals <- creosote.summary %>%
    group_by(source) %>%
    summarize(med = median(sum_chi2)) %>%
    arrange(med)
  creosote.summary$source <- factor(creosote.summary$source, levels = order.vals$source)
  
  p <- ggplot(creosote.summary, aes(x = source, y = sum_chi2)) +
    geom_boxplot() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = '', y = 'Sum Chi2')
  
  ggsave('6_profile_data/doc/chi2_allsites_withcreosote.png', p)
  
    
  
  
  
  
}