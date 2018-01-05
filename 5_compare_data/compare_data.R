merge_studies <- function() {
  mke <- raw_mke
  # first, handle zeros in mke data by turning values below dl to 0
  mke$remark_cd[is.na(mke$remark_cd)] <- "d"  # if remark_cd is blank it's a detection
  mke$remark_cd <- ifelse(mke$remark_cd == "E", "d", mke$remark_cd)  # classify Estimated values as detections
  mke$valuetouse <- ifelse(mke$remark_cd == "d", mke$result_va, 0)
  mke.m <- select(mke, site_no, parm_cd, valuetouse) %>%
    rename(pcode = parm_cd)
  
  # merge MKE samples with GLRI samples
  glri <- samples
  glri.m <- mutate(glri, site_no = paste0("0", STAID)) %>%
    select(site_no, pcode, RESULT, PARAM_SYNONYM, EPApriority16)
  
  ###########################
  # get pilot study from MKE
  pilot <- read_excel('M:/QW Monitoring Team/GLRI toxics/GLRI II/WY2017/WY 2017 planning/PilotTest2016/Data/Sediment data/GLRI2016PilotSediment.xlsx',
                      skip = 8, sheet = "Main")
  pilot.m <- pilot[-(1:14), -2]
  pilot.m <- pilot.m[-(39:47),]
  
  grep("surrogate", pilot.m$`Client ID`, ignore.case = T)
  names(pilot.m)[seq(2, 42, 2)] <- gsub("-MIR-SED", "", names(pilot.m)[seq(2, 42, 2)])
  names(pilot.m)[seq(3, 43, by = 2)] <- paste0(names(pilot.m)[seq(2, 42, 2)], "-R")
  vals <- select(pilot.m, 1, seq(2, 42, 2)) %>%
    gather(site_id, value, -`Client ID`)
  
  remarks <- select(pilot.m, 1, seq(3, 43, by = 2))
  names(remarks)[2:22] <- names(pilot.m)[seq(2, 42, 2)]
  remarks <- gather(remarks, site_id, remark, -`Client ID`)
  names(vals)[1] <- "compound"
  names(remarks)[1] <- "compound"
  
  pilot.merged <- left_join(vals, remarks)
  pilot.merged <- mutate(pilot.merged, STAT_ID = gsub("-A|-B|-C", "",site_id))
  
  #find MKE sites from glri dataset and merge to get site numbers
  mke.sites <- glri[grep("milwaukee", glri$Watershed, ignore.case = T), ]
  mke.sites <- select(mke.sites, STAT_ID, STAID) %>%
    distinct() %>%
    rename(site_no = STAID)
  
  # get pcodes and compound names from glri to merge with pilot
  compound.m <- select(glri, PARAM_SYNONYM, pcode) %>% 
    distinct() %>%
    rename(compound = PARAM_SYNONYM)
  
  pilot.merged <- left_join(pilot.merged, mke.sites)
  pilot.merged <- left_join(pilot.merged, compound.m)
  
  # merge datasets based on site number and pcode
  df <- inner_join(mke.m, glri.m, by = c('site_no', 'pcode')) %>%
    rename("GLRI_conc_ppb" = "RESULT",
           "MKE_conc_ppb" = "valuetouse",
           "Compound" = "PARAM_SYNONYM")
  
  
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