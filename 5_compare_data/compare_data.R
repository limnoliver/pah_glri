library(dplyr)
library(tidyr)
library(ggplot2)
library(dataRetrieval)
merge_studies <- function() {
  #############################
  # get MKE study
  mke.m <- select(processed_mke, site_no, parm_cd, MKE) %>%
    rename(pcode = parm_cd)
  
  ##############################
  # get GLRI samples
  
  glri <- samples
  glri.m <- mutate(glri, site_no = ifelse(nchar(STAID) != 15, paste0("0", STAID), STAID)) %>%
    select(site_no, pcode, RESULT, PARAM_SYNONYM) %>%
    rename(GLRI = RESULT, study_compound_name = PARAM_SYNONYM)
  
  ###########################
  # get PILOT study data
  pilot <- read_excel('M:/QW Monitoring Team/GLRI toxics/GLRI II/WY2017/WY 2017 planning/PilotTest2016/Data/Sediment data/GLRI2016PilotSediment.xlsx',
                      skip = 8, sheet = "Main")
  pilot <- pilot[-(1:14), -2]
  pilot.m <- pilot[-(39:47),]
  
  names(pilot.m)[seq(2, 42, 2)] <- gsub("-MIR-SED", "", names(pilot.m)[seq(2, 42, 2)])
  names(pilot.m)[seq(3, 43, by = 2)] <- paste0(names(pilot.m)[seq(2, 42, 2)], "-R")
  vals <- select(pilot.m, 1, seq(2, 42, 2)) %>%
    gather(site_id, value, -`Client ID`)
  vals$value <- as.numeric(vals$value)
  
  remarks <- select(pilot.m, 1, seq(3, 43, by = 2))
  names(remarks)[2:22] <- names(pilot.m)[seq(2, 42, 2)]
  remarks <- gather(remarks, site_id, remark, -`Client ID`)
  names(vals)[1] <- "compound"
  names(remarks)[1] <- "compound"
  
  pilot.m <- left_join(vals, remarks)
  pilot.m <- mutate(pilot.m, STAT_ID = gsub("-A|-B|-C", "",site_id)) %>%
    rename(PILOT = value)
  
  # if remark = U, the value is censored
  pilot.m$PILOT <- ifelse(pilot.m$remark %in% "U", 0, pilot.m$PILOT)
  pilot.m <- select(pilot.m, -remark)
  
  #find MKE sites from glri dataset and merge to get site numbers
  mke.sites <- glri[grep("milwaukee", glri$Watershed, ignore.case = T), ]
  mke.sites <- select(mke.sites, STAT_ID, STAID) %>%
    distinct() %>%
    rename(site_no = STAID)
  
  # get pcodes and compound names from glri to merge with pilot
  compound.m <- select(glri, PARAM_SYNONYM, pcode) %>% 
    distinct() %>%
    rename(compound = PARAM_SYNONYM)
  
  pilot.m <- left_join(pilot.m, mke.sites)
  pilot.m <- left_join(pilot.m, compound.m)
  
  # because there are triplicates, and we'll want to compare individual samples
  # and the mean of samples, create a mean value that will have the original site_no,
  # but then the triplicate samples have an a-b-c attached to the site no so 
  # there is no attempt to double merge those samples
  
  pilot.avg <- group_by(pilot.m, site_no, pcode, compound) %>%
    summarize(PILOT = mean(PILOT))
  
  pilot.m <- mutate(pilot.m, 
                    site_no = paste0(site_no, gsub("[[:alpha:]]{3}-", "", site_id)))
  
  pilot.m <- bind_rows(pilot.m, pilot.avg) %>%
    select(site_no, pcode, PILOT)
  
  pilot.m <- mutate(pilot.m, site_no = ifelse(nchar(site_no) != 15, paste0("0", site_no), site_no))
  
  
  # get GLRI data that measured six PAH compounds at all sites
  # as part of schedule 5433
  
  # merge datasets based on site number and pcode
  # should keep all possible site numbers to be filtered out 
  # later once all data are merged
  df <- full_join(mke.m, glri.m, by = c('site_no', 'pcode'))
  
  df <- full_join(df, pilot.m, by = c('site_no', 'pcode'))
  
  # now merge in site/compound metadata 
  names(glri)
  site.dat <- select(glri, "State", "Watershed", "Lake", "STAID", "Site", "STAT_ID", "Lat", "Lon")
  site.dat <- mutate(site.dat, 
                     site_no = ifelse(nchar(STAID) != 15, paste0("0", STAID), STAID)) %>%
    select(-STAID) %>%
    rename(site_code = STAT_ID)
  
  # get as much site info as possible by merging with glri dataset which has most sites
  df <- left_join(df, site.dat, by = 'site_no') 
  
  # now merge in compound metadata
  compound.dat <- select(glri, )
  # write in merge of glri6 data when it becomes available
  return(df)
  
  
  p <- ggplot(df, aes(x = GLRI_conc_ppb, y = MKE_conc_ppb)) +
    geom_point(aes(color = site_no), alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0)  + 
    labs(x = "GLRI Concentration (ppb)", y = "MKE Study Concentration (ppb)") +
    theme_bw() +
    scale_x_continuous(breaks = c(0, 10, 100, 1000, 10000), trans = "log1p") +
    scale_y_continuous(breaks = c(0, 10, 100, 1000, 10000), trans = "log1p") +
    geom_text(label = "1:1 Line",aes(x = 18, y = 8), color = "black")
  ggsave("GLRI_MKE_comparison_scatter.png", p)
  
  df.long <- gather(df, study, value, -site_no, -pcode, -Compound, -EPApriority16)
  df.long <- mutate(df.long, EPA = ifelse(EPApriority16 == TRUE, "EPA Priority", "Other PAHs"))
  p <- ggplot(df.long, aes(x = Compound, y = value)) +
    geom_point(aes(color = study), alpha = 0.8) +
    facet_grid(site_no~EPA, scale = "free_x", margins = F) +
    scale_y_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 10000)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Concentration (ppb)")
  ggsave("GLRI_MKE_comparison_bysite.png", p)
  
  # calculate percent difference in each site/compound
  head(df)
  df.diff <- mutate(df, normalized_diff = (MKE_conc_ppb - GLRI_conc_ppb)/mean(c(MKE_conc_ppb, GLRI_conc_ppb))*100) %>%
    mutate(EPA = ifelse(EPApriority16 == TRUE, "EPA Priority", "Other PAHs"))
  
  p <- ggplot(df.diff, aes(x = Compound, y = normalized_diff)) +
    geom_boxplot() +
    geom_hline(yintercept = 0, color = "red", alpha = 0.8) +
    geom_text(label = "MKE > GLRI", data = data.frame(Compound = "Acenaphthylene", 
                                                      normalized_diff = 300,
                                                      EPA = "EPA Priority"),
              color = "red") +
    geom_text(label = "MKE < GLRI", data = data.frame(Compound = "Acenaphthylene", 
                                                      normalized_diff = -100,
                                                      EPA = "EPA Priority"),
              color = "red") +
    facet_grid(~EPA, scale = "free_x") +
    labs(y = "Normalized Difference Between Studies \n(Milwaukee - GLRI/mean)", x = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("GLRI_MKE_diff_bycompound.png", p)
  
  p <- ggplot(df.diff, aes(x = site_no, y = normalized_diff)) +
    geom_boxplot(aes(fill = EPA)) +
    #facet_grid(EPA~., scale = "free_x") +
    geom_text(label = "MKE > GLRI", data = data.frame(site_no = "04087040", 
                                                      normalized_diff = 500),
              color = "red") +
    geom_text(label = "MKE < GLRI", data = data.frame(site_no = "04087040", 
                                                      normalized_diff = -100),
              color = "red") +
    geom_hline(yintercept = 0, color = "red", alpha = 0.8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "Normalized Difference Between Studies \n(Milwaukee - GLRI/mean)",
         x = "Site")
  
  ggsave("GLRI_MKE_diff_bysite.png", p)
  
}
library(readxl)
merge_recovery <- function() {

  # get recovery data from MKE (NWQL) and GLRI (Batelle)
  mke_recovery <- read.csv('5_compare_data/raw/NWQL_prepSpike_recoveries.csv', 
                           skip = 2, stringsAsFactors = F)
  head(mke_recovery)
  
  # get glri recover data
  file.loc <- "M:/QW Monitoring Team/GLRI toxics/GLRI II/WY2017/Data"
  sample.files <- list.files("M:/QW Monitoring Team/GLRI toxics/GLRI II/WY2017/Data")
  sample.files <- grep("S17", sample.files, value = T)
  sample.files <- grep(".xlsx", sample.files, value = T)
  
  all_dat <- data.frame()
  for (i in sample.files) {
    temp_dat <- read_excel(file.path(file.loc, i), sheet = "MS", 
                           skip = 21, col_names = T)
    temp_dat_filt <- temp_dat[,c(1, grep('recovery', names(temp_dat), ignore.case = T), grep('qualifier', names(temp_dat), ignore.case = T))]
    temp_dat_filt <- temp_dat_filt[grep("^Naphthalene$", temp_dat_filt[[1]]):grep("perylene", temp_dat_filt[[1]]),]
    # get rid of estimates in samples where the matrix spike was 
    # < 5x greater than the sample conc
    temp_dat_filt$`% Recovery` <- ifelse(temp_dat_filt$Qualifier %in% "n", NA, temp_dat_filt$`% Recovery`)
    
    all_dat <- bind_cols(temp_dat_filt, all_dat)
    
  }
  # now find median, min, max
  all_dat <- all_dat[,c(1,2,5,8,11)]
  all_dat$min <- apply(all_dat[,2:5], 1, FUN = min, na.rm = T)
  all_dat$min[all_dat$min == Inf] <- NA
  all_dat$mean <- apply(all_dat[,2:5], 1, FUN = mean, na.rm = T)
  all_dat$mean[is.nan(all_dat$mean)] <- NA
  all_dat$max <- apply(all_dat[,2:5], 1, FUN = max, na.rm = T)
  all_dat$max[all_dat$max == -Inf] <- NA
  all_dat$type <- "spikes"
  all_dat$source <- "Battelle"
  
  glri_spikes <- all_dat %>%
    select(Units, min, mean, max, type, source) %>%
    rename(compound = Units) %>%
    mutate(compound_simple = tolower(gsub("[[:punct:]]", "", compound))) %>%
    filter(!is.na(mean))
  
  mke_spikes <- select(mke_recovery, Parameter, Minimum, Mean, Maximum) %>%
    rename(compound = Parameter, min = Minimum, mean = Mean, max = Maximum) %>%
    mutate(compound_simple = tolower(gsub("[[:punct:]]", "", compound))) %>%
    filter(compound_simple %in% glri_spikes$compound_simple)
  
  mke_spikes$type <- 'spikes'
  mke_spikes$source <- 'NWQL'


  spikes <- bind_rows(mke_spikes, glri_spikes) %>%
    filter(!is.na(mean)) %>%
    filter(compound_simple != 'biphenyl')
  
  compound.order <- filter(spikes, source == "Battelle") %>%
    arrange(mean) %>%
    select(compound_simple)
  spikes$compound_simple <- factor(spikes$compound_simple, levels = compound.order[[1]])
  
  
  # plot data
  p <- ggplot(spikes, aes(x =compound_simple, y = mean, group = source)) +
    geom_point(aes(color = source), size = 4, shape = 15, position = position_dodge(width =0.6)) +
    geom_pointrange(aes(ymin = min, ymax = max, color = source), position = position_dodge(width =0.6)) +
    theme_bw() +
    theme(panel.grid.minor.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(y = "% Recovery from Spikes", x = "", 
         title = "Comparison of matrix spike recoveries from NWQL and Battelle",
         subtitle = "Values represent the minimum, mean, and maximum reported recoveries. 
Maximum n value for any compound Battelle = 4, NWQL = 26.")

  ggsave('spike_comparison.png', p)
  
  #########################################################
  # this glri data is actually surrogates - find mke surrogates
  glri_sur_m <- samples %>%
    filter(UNIT == "PCT_REC") %>%
    select(PARAM_SYNONYM, RESULT) %>%
    group_by(PARAM_SYNONYM) %>%
    summarize_at(vars(RESULT), funs(min, median, max), na.rm = T)
  
  glri_sur_q <- samples %>%
    filter(UNIT == "PCT_REC") %>%
    select(PARAM_SYNONYM, RESULT) %>%
    group_by(PARAM_SYNONYM) %>%
    summarize(first = quantile(RESULT, probs = 0.25),
              third = quantile(RESULT, probs = 0.75))
  
  glri_sur <- left_join(glri_sur_m, glri_sur_q) 
  glri_sur_names <- glri_sur$PARAM_SYNONYM
  
  glri_sur <- select(glri_sur, min, first, median, third, max) %>%
    t()
  
  compounds <- gsub("-[[:alnum:]]+", "", glri_recovery$PARAM_SYNONYM)
  compounds <- gsub("\\(", "\\[", compounds)
  compounds <- gsub("\\)", "\\]", compounds)

  mke_sur <- processed_mke %>% 
    rename(parameter_cd = parm_cd) %>%
    left_join(parameterCdFile)
  row.keep <- grep("recovery", mke_sur$parameter_nm)
  mke_sur <- mke_sur[row.keep, ] 
  
  mke_sur_m <- group_by(mke_sur, parameter_cd) %>%
    summarize_at(vars(result_va), funs(min, median, max), na.rm = T)
  
    
  
  recovery <- matrix(ncol = 8, nrow = 5)
  recovery[,1:8] <- c(mke_rec[,1], glri_rec[,1],
                      mke_rec[,2], glri_rec[,2],
                      mke_rec[,3], glri_rec[,3],
                      mke_rec[,4], glri_rec[,4]) 
  test <- data.frame(`Percent Recovery` = c(1:24), 
                     Lab = rep(c("A", "B"), 12), 
                     Chemicals = c(rep("chem1", 6), rep("chem2", 6),
                                   rep("chem3", 6), rep("chem4", 6)))
  test.p <- boxplot(test$Percent.Recovery ~ test$Lab*test$Chemicals, col = 'red')
  test.p$stats <- recovery
  test.p$names <- c(compounds[1], "", compounds[2], "", compounds[3], "", compounds[4], "")
  
  png("Method_pctrecovery_comparison.png", height = 500, width = 1000)
  par(mar=c(5,5,2,2))
  bxp(z = test.p, border = c("black", "red"), outcol = "black", 
      ylab = "Percent Recovery", ylim = c(40, 145), xaxt = 'n', cex.axis = 2, cex.lab = 2)
  legend('topleft', legend = c("NWQL", "Batelle"), border = c('black', 'red'), 
         cex = 2, fill = "white")
  text(x = c(1.5, 3.5, 5.5, 7.5), y = 25, labels = compounds, xpd = T, cex = 2)
  dev.off()
  # create a boxplot out of the following values: min, 1st quartile, median, 2nd quartile, max
}
