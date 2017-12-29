# script to process original raw data
# this is based on a number of notes made by Pete Lenaker in 
# the file: 'M:/QW Monitoring Team/GLRI toxics/GLRI II/WY2017/Data/GLRI PAH Battelle Data Notes_READ ME'

process_sample_merge <- function(raw_sample, raw_site){
  dat.clean <- raw_sample
  dat.clean$STAT_ID[dat.clean$STAT_ID == "UJC"] <- "UCJ" # fix wrong site id for Underwood Creek
  dat.clean$STAT_ID[dat.clean$STAT_ID == "CRO"] <- "CRP" # fix wrong site id for Underwood Creek
  
  # create just a date (collection) column without time
  dat.clean$COLLECTION_DATE_NOTIME <- as.Date(dat.clean$COLLECTION_DATE)
  
  # create two datasets, one where there should be no issue merging by site,
  # and the other where STAT_ID is duplicated across states
  dat.clean.nodup <- filter(dat.clean, !(STAT_ID %in% c("CRM", "RRS", "RRR", "RRB")))
  dat.clean.dup <- filter(dat.clean, (STAT_ID %in% c("CRM", "RRS", "RRR", "RRB")))
  
  dat.clean.nodup <- left_join(dat.clean.nodup, raw_site, by = c("STAT_ID" = "Base.Field.ID"))
  
  # fix spots where a STAT_ID was used by multiple states
  #CRM
  dat.clean.dup$State <- NA
  dat.clean.dup$State[dat.clean.dup$STAT_ID == 'CRM' & 
                        dat.clean.dup$COLLECTION_DATE_NOTIME == as.Date("2017-06-28")] <- "MI"
  dat.clean.dup$State[dat.clean.dup$STAT_ID == 'CRM' & 
                        dat.clean.dup$COLLECTION_DATE_NOTIME == as.Date("2017-06-07")] <- "OH"
  #RRS
  dat.clean.dup$State[dat.clean.dup$STAT_ID == 'RRS' & 
                        dat.clean.dup$COLLECTION_DATE_NOTIME == as.Date("2017-06-28")] <- "MI"
  dat.clean.dup$State[dat.clean.dup$STAT_ID == 'RRS' & 
                        dat.clean.dup$COLLECTION_DATE_NOTIME == as.Date("2017-06-06")] <- "OH"
  #RRR
  dat.clean.dup$State[dat.clean.dup$STAT_ID == 'RRR' & 
                        dat.clean.dup$COLLECTION_DATE_NOTIME == as.Date("2017-06-28")] <- "MI"
  dat.clean.dup$State[dat.clean.dup$STAT_ID == 'RRR' & 
                        dat.clean.dup$COLLECTION_DATE_NOTIME == as.Date("2017-06-07")] <- "WI"
  #RRB
  dat.clean.dup$State[dat.clean.dup$STAT_ID == 'RRB' & 
                        dat.clean.dup$COLLECTION_DATE_NOTIME == as.Date("2017-06-28")] <- "MI"
  dat.clean.dup$State[dat.clean.dup$STAT_ID == 'RRB' & 
                        dat.clean.dup$COLLECTION_DATE_NOTIME == as.Date("2017-06-06")] <- "OH"
  
  # now merge the cleaned up duplicates by STATID and state
  dat.clean.dup <- left_join(dat.clean.dup, raw_site, by =c("STAT_ID" = "Base.Field.ID", "State"))
  
  # append dat.clean.dup to dat.clean.nodup
  dat.clean <- bind_rows(dat.clean.nodup, dat.clean.dup)
  
  # get rid of samples that are "bad"
  dat.clean <- filter(dat.clean, LAB_SAMPLE_ID != "K7216") # get rid of first sample taken at SLR
  
  # get rid of columns that are unnecessary
  dat.clean <- select(dat.clean, -c(1:7, 14:21))
  
  # merge with compound metadata
  dat.clean2 <- get_compound_info(pah_dat = dat.clean, merge_type = 'name', merge_col = 'PARAM_SYNONYM')
  
  # export dat.clean
  return(dat.clean2)
  
}

process_nondetects <- function(processed_sample_merge, detect.decision = "zero",
                               code.column = "LAB_QUAL", detect.code = "U", 
                               detect.limit.column = "DETECT_LIMIT") {
  if (detect.decision == 'zero') {
    processed_sample <- 
      mutate(processed_sample_merge, 
             RESULT = ifelse(grepl(detect.code, processed_sample_merge[[code.column]]), 0, RESULT))
  } else if (detect.decision == 'half'){
    processed_sample <- 
      mutate(processed_sample_merge, 
             RESULT = ifelse(grepl(detect.code, processed_sample_merge[[code.column]]), 0.5*processed_sample_merge[[detect.limit.column]], RESULT))
  }
 
}

