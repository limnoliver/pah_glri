# return blanks and corresponding samples with same state/site id/date
# this will return samples, duplicates, and blanks for sites with blanks
filter_blanks <- function(processed_sample) {
  # create a unique id for state/site/date combination
  processed_sample$unique_id <- paste(processed_sample$State, processed_sample$STAT_ID, processed_sample$COLLECTION_DATE_NOTIME, sep = "_")
  
  # find blanks and get the unique id of all blanks
  blank_sample <- filter(processed_sample, processed_sample$FIELD_QC_CODE == "FB")
  blank_sample_ids <- unique(blank_sample$unique_id)
  
  # get all samples from the unique IDS of blanks
  # this will also return duplicates and samples from the same state/site/date to compare the blanks to
  corresponding_samples <- filter(processed_sample, unique_id %in% blank_sample_ids)
  return(corresponding_samples)
}

# return duplicates
# this will return duplicates and corresponding samples from the same state/site/date

filter_duplicates <- function(processed_sample) {
  processed_sample$unique_id <- paste(processed_sample$State, processed_sample$STAT_ID, processed_sample$COLLECTION_DATE_NOTIME, sep = "_")
  
  # find duplicates and get the unique id of all duplicates
  duplicate_sample <- filter(processed_sample, processed_sample$FIELD_QC_CODE == "DU")
  duplicate_sample_ids <- unique(duplicate_sample$unique_id)
  corresponding_samples <- filter(processed_sample, unique_id %in% duplicate_sample_ids) %>%
    filter(FIELD_QC_CODE != "FB")
  return(corresponding_samples)
}

# filter out duplicates and blanks

filter_samples <- function(processed_sample) {
  samples <- filter(processed_sample, FIELD_QC_CODE == "SA")
  return(samples)
}

filter_surrogates <- function(processed_sample) {
  glri_sur_m <- processed_sample %>%
    filter(UNIT == "PCT_REC") %>%
    select(PARAM_SYNONYM, RESULT)
  return(glri_sur_m)
}

filter_surrogates_5507 <- function(processed_5507) {
  
  mke <- processed_5507 %>% 
    rename(parameter_cd = parm_cd) %>%
    left_join(parameterCdFile)
  
  row.keep <- grep("recovery", mke$parameter_nm)
  mke <- mke[row.keep, ] 
  
  return(mke)
  
}

filter_surrogates_5433 <- function(processed_5433) {
  
  s_5433 <- processed_5433 %>% 
    rename(parameter_cd = parm_cd) %>%
    left_join(parameterCdFile)
  
  row.keep <- grep("recovery", s_5433$parameter_nm)
  s_5433 <- s_5433[row.keep, ] 
  
  return(s_5433)
  
}

filter_characteristics <- function(samples) {
  sed_char <- filter(samples, CLASS %in% c('SOLIDS', 'TOC', 'GS_P')) %>%
    select(unique_site_id = unique_id, state = State, watershed = Watershed, lake = Lake, 
           usgs_station_id = STAID, site_description = Site, latitude = Lat, longitude = Lon, 
           collection_date = COLLECTION_DATE_NOTIME, parameter = PARAM_SYNONYM, 
           parameter_description = DESCR, result = RESULT, lab_qualifier = LAB_QUAL, 
           result_unit = UNIT, detection_limit = DETECT_LIMIT, detection_limit_code = DETECT_LIMIT_CODE)
}
