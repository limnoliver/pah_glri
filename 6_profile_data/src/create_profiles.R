# calculate difference between source and sample profiles
library(dplyr)
library(rlang)
library(tidyr)
library(pah)

# remove all samples where there is a zero - meaning BDL
prep_profiles <- function(pah_dat) {
  pah_dat$sample_id <- paste0(pah_dat$State, "-", pah_dat$STAT_ID)
  pah_dat <- filter(pah_dat, sourceProfile12 == TRUE) # keep only chemicals in source profiles
  pah_dat.remove <- unique(pah_dat$sample_id[pah_dat$RESULT == 0]) # get rid of pah_dat BDL
  pah_dat <- filter(pah_dat, !(sample_id %in% pah_dat.remove))
  pah_dat <- select(pah_dat, sample_id, RESULT, casrn, Lake, State, Watershed) #narrow columns
  return(pah_dat)
}

# function to pull out pah 16 compounds and order samples by totals
# then group totals into bins for later plotting
prep_totals <- function(pah_dat) {
  pah_dat$sample_id <- paste0(pah_dat$State, "-", pah_dat$STAT_ID)
  samples_16 <- filter(pah_dat, PARAM_SYNONYM == 'Priority Pollutant PAH') %>%
    arrange(RESULT)
  samples_16$sample_id <- paste0(samples_16$State, "-", samples_16$STAT_ID)
  samples_16 <- select(samples_16, sample_id, RESULT) %>%
    group_by(sample_id) %>%
    summarize(Priority16 = mean(RESULT)) %>%
    mutate(Priority16_bin = ntile(Priority16, 4))
  
  # summarize by bin
  summary_16 <- group_by(samples_16, Priority16_bin) %>%
    summarize(round(max(Priority16), 0))
  
  samples_16$Priority16_bin[samples_16$Priority16_bin == 1] <- "Low (< 443 ppb)"
  samples_16$Priority16_bin[samples_16$Priority16_bin == 2] <- "Med (443-2489 ppb)"
  samples_16$Priority16_bin[samples_16$Priority16_bin == 3] <- "Med (2489-13426 ppb)"
  samples_16$Priority16_bin[samples_16$Priority16_bin == 4] <- "High (>13426 ppb)"
  
  # order the factor
  samples_16$Priority16_bin <- factor(samples_16$Priority16_bin, levels = c("Low (< 443 ppb)","Med (443-2489 ppb)",
                                                                            "Med (2489-13426 ppb)","High (>13426 ppb)"))
  return(samples_16)
}

