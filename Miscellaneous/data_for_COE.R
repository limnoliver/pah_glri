# grab pah and owc data for COE collaborators
library(dplyr)
pah <- make('samples') %>%
  select(date = COLLECTION_DATE_NOTIME, station_id = STAID, site_abbrev = SITE, 
         site_description = Site, state = State, watershed = Watershed, 
         analysis_method = ANALYSIS_METH, cas_num = CAS_NUM, parameter = PARAM_SYNONYM, parameter_description = DESCR, 
         parameter_class = CLASS, result = RESULT, lab_qualifier = LAB_QUAL, unit = UNIT,
         detect_limit = DETECT_LIMIT, detect_limit_code = DETECT_LIMIT_CODE)

p5433 <- make('processed_5433') %>%
  left_join(dataRetrieval::parameterCdFile, by = c('parm_cd' = 'parameter_cd')) %>%
  select(date = sample_dt, station_id = site_no, analysis_method = meth_cd,
         cas_num = casrn, parameter = srsname, parameter_description = parameter_nm, 
         parameter_class = parameter_group_nm, result = conc_5433, lab_qualifier = remark_cd,
         units = parameter_units)

write.csv(pah, file = 'Miscellaneous/pah_GLRI_2017_COE.csv', row.names = F)
write.csv(p5433, file = 'Miscellaneous/owc_GLRI_2017_COE.csv', row.names = F)

  