# Analysis of 12-compound PAH profiles

library(remake)
library(dplyr)
library(pah)
library(Hmisc)
library(reshape2)


profiles <- make("profiles")

chi <- profiles$sum_chi2

pro <- profiles$profiles

# Sum chi2 values for profiles
profilesChi2 <- profiles$sum_chi2

# summarize chi2 values for each source
chiSummary <- group_by(chi, source)%>%
  summarise(median= median(sum_chi2, na.rm = TRUE))

# output summary to csv
write.csv(chiSummary, "C:/Users/akbaldwi/Documents/GLRI/PAHs/pah/pah_glri/austin_analysis/profiles_sumChi2_medianBySource.csv", row.names = F)


#==============================================================
# how similar are profiles from site to site?
#   create correlation matrix to compare sites
#   rcorr method (Hmisc package)

proCorr <- select(pro, unique_id, abbrev=Abbreviation, prop_conc)

# get unique rows 
proCorr <- unique(proCorr)

proCorr2 <- dcast(proCorr, abbrev ~ unique_id, value.var="prop_conc")

# Correlation matrix: 
corr <- select(proCorr2, 2:71)
corr <- as.matrix(corr)

# Spearman:
profileCorrSpearman <- rcorr(corr, type="spearman")

# correlation coefficient
profileCorrSpearman_Corr <- as.data.frame(profileCorrSpearman[1])
profileCorrSpearman_Corr <- cbind(compound = rownames(profileCorrSpearman_Corr), profileCorrSpearman_Corr)
rownames(profileCorrSpearman_Corr) <- NULL

# p-value
profileCorrSpearman_Pval <- as.data.frame(profileCorrSpearman[3])
profileCorrSpearman_Pval <- cbind(compound = rownames(profileCorrSpearman_Pval), profileCorrSpearman_Pval)
rownames(profileCorrSpearman_Pval) <- NULL

# save correlation coefficients and pvalues
write.csv(profileCorrSpearman_Corr, "C:/Users/akbaldwi/Documents/GLRI/PAHs/pah/pah_glri/austin_analysis/Spearman_correlations_between_12compound_profiles_across_sites.csv", row.names = F)
write.csv(profileCorrSpearman_Pval, "C:/Users/akbaldwi/Documents/GLRI/PAHs/pah/pah_glri/austin_analysis/Spearman_pValues_between_12compound_profiles_across_sites.csv", row.names = F)

#=======================================================================
# plot profiles of samples together to show similarity

# output profiles to id those with unusual profiles
write.csv(proCorr2, "C:/Users/akbaldwi/Documents/GLRI/PAHs/pah/pah_glri/austin_analysis/profiles_for_each_site.csv", row.names = F)

setwd("C:/Users/akbaldwi/Documents/GLRI/PAHs/pah/pah_glri/austin_analysis/")

# import table with parent/alkyl and HMW/LMW ratios.
#   use ratios to differentiate profiles on plot (petrogenic vs pyrogenic by color)
ratios <- read.csv("C:/Users/akbaldwi/Documents/GLRI/PAHs/pah/pah_glri/9_parent_weight/doc/summary_table.csv")

# create pyrogenic vs petrogenic column, based on ratios
ratios$pyroPetro <- ifelse(ratios$parent_alkyl <1 & ratios$HMW_LMW <1, "petrogenic", "pyrogenic")
ratios <- select(ratios, unique_id, pyroPetro)

# merge pyrogenic/petrogenic column into profiles df
proCorr <- merge(proCorr, ratios, by="unique_id")

# set order of factors for plotting (low to high mol wt)
proCorr$abbrev <- factor(proCorr$abbrev, levels=c("Phen","Anth","FluA","Pyr","BaA","Ch","BbF","BkF","BeP","BaP","IndPy","BghiP"))

ggplot(proCorr, aes(x=abbrev, y=prop_conc, group=unique_id, color=pyroPetro))+
  geom_line(data=filter(proCorr, pyroPetro=="pyrogenic"),alpha=0.75)+
  geom_point(data=filter(proCorr, pyroPetro=="pyrogenic"),size=1, alpha=0.7)+
  geom_line(data=filter(proCorr, pyroPetro=="petrogenic"),alpha=0.85)+
  geom_point(data=filter(proCorr, pyroPetro=="petrogenic"),shape=21,fill="white",size=1,alpha=0.70)+
  xlab("Individual PAH compounds\n(ordered from low to high molecular wt)") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position=c(0.9,0.9), legend.justification = c(0.9,0.9),
        legend.title=element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=rel(.8), angle=90, hjust=1, vjust=0.5),
        axis.text.y = element_text(size=rel(.8)))
ggsave("profiles__of_sites.pdf", width=3.5, height=3., dpi=300)
ggsave("profiles__of_sites.png", width=3.5, height=3., dpi=300)
shell.exec("profiles__of_sites.pdf")







