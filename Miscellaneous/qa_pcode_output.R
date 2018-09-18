# pcodes for Pete
samples <- make('samples')

library(dplyr)
compounds <- select(samples, parameter = PARAM_SYNONYM, CAS_NUM, pcode, battelle_method = ANALYSIS_METH, battelle_unit =UNIT) %>%
  distinct() %>%
  slice(17:53)

compounds$pcode[compounds$parameter == "Biphenyl"] <- '63752'
  
compounds <- left_join(compounds, dataRetrieval::parameterCdFile, by = c('pcode' = 'parameter_cd'))

write.csv(compounds, 'Miscellaneous/compounds_with_pcodes.csv', row.names = F)

sites <- distinct(select(samples, unique_id)) 

wi_bro <- filter(samples, unique_id == "WI-BRO")
wi_bro$COLLECTION_DATE
mn_slr <- filter(samples, unique_id == 'MN-SLR')

# write QA data for Pete
lab_duplicates <- make('lab_duplicates')
write.csv(lab_duplicates, 'Miscellaneous/battelle_lab_duplicates.csv', row.names = F)

proc_blanks <- make('procedural_blanks')
write.csv(proc_blanks, 'Miscellaneous/battelle_procedural_blanks.csv', row.names = F)

pct_rec_mspikes <- make('pct_rec_mspikes')
write.csv(pct_rec_mspikes, 'Miscellaneous/battelle_matrixspikes_pct_recovery.csv', row.names = F)

pct_rec_labcontrol <- make('pct_rec_labcontrol')
write.csv(pct_rec_labcontrol, 'Miscellaneous/battelle_labcontrol_pct_recovery.csv', row.names = F)
