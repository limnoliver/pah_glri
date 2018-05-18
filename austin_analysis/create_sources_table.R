# create a table similar to that in Baldwin et al 2017 (table 1)
# that shows all sources, including source category, source description with abbreviation, 
# which analysis methods that source was considered in, and referenes

# the "methods" column should include DR = diagnostic ratios, 
# P = profiles, PCA, and CMB 

# so the DR should come from the table pah::source_ratios
# the P should come from pah::source_profiles
# and merged with info from pah::sources

# get the sources (by abbrev) that were considered in the profiles
profiles <- make('profiles')
sources_in_profs <- unique(profiles$sum_chi2$source[!is.na(profiles$sum_chi2$sum_chi2)])
sources_in_profs <- data.frame(source_abbreviation = sources_in_profs,
                               method_p = 'P')

# get the sources that were used in diagnostic ratios
ratios <- make('ratio_distance')
sources_in_ratios <- ratios$source$source
sources_in_ratios <- data.frame(source_abbreviation = sources_in_ratios,
                                method_dr = 'DR')

# source_table
sources_table <- full_join(sources_in_ratios, sources_in_profs)
sources_table <- left_join(sources_table, pah::source_ratios[, c('abbrev', 'Reference')],
                           by = c('source_abbreviation' = 'abbrev'))
sources_table$methods <- paste(sources_table$method_dr, sources_table$method_p, sep = ', ')
sources_table$methods <- gsub('NA,\\s|,\\sNA', '', sources_table$methods)

# now merge with source info to get sources category

table <- left_join(sources_table, pah::sources[, c('source_subcategory', 'source_abbrev', 'source_short_no_ref')], 
                                                       by = c('source_abbreviation' = 'source_abbrev')) %>%
  mutate(name = paste0(source_short_no_ref, ' (', source_abbreviation, ')')) %>%
  select(source_subcategory, name, Reference, methods)

