library(dplyr)
library(tidyr)
merge_studies <- function() {
  #############################
  # get MKE study
  mke <- raw_mke
  # first, handle zeros in mke data by turning values below dl to 0
  mke$remark_cd[is.na(mke$remark_cd)] <- "d"  # if remark_cd is blank it's a detection
  mke$remark_cd <- ifelse(mke$remark_cd == "E", "d", mke$remark_cd)  # classify Estimated values as detections
  mke$MKE<- ifelse(mke$remark_cd == "d", mke$result_va, 0)
  mke.m <- select(mke, site_no, parm_cd,MKE) %>%
    rename(pcode = parm_cd)
  
  ##############################
  # get GLRI samples
  
  glri <- samples
  glri.m <- mutate(glri, site_no = ifelse(nchar(STAID) != 15, paste0("0", STAID), STAID)) %>%
    select(site_no, pcode, RESULT) %>%
    rename(GLRI = RESULT)
  
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