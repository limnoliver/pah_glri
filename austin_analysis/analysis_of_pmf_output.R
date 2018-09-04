# Script to compare PMF output profiles to source profiles using plots and chi-squared.
#=========================================================================================

pathToSave <- paste("C:/Users/akbaldwi/Documents/EPA PMF/Output/", sep="")

library(reshape2)
library(ggplot2)
library(dplyr)
library(matrixStats)
library(RColorBrewer)
#====================================================================

# import PMF factors
pmf <- read.csv("C:/Users/akbaldwi/Documents/EPA PMF/Output/glri11_2factor_summary_for_R.csv", stringsAsFactors=FALSE)
# OR
pmf <- read.csv("C:/Users/akbaldwi/Documents/EPA PMF/Output/glri10_3factor_summary_for_R.csv", stringsAsFactors=FALSE)

# import source profiles from literature
sources <- read.csv("C:/Users/akbaldwi/Documents/MMSD_PAHs/R/for_GLRI_reuse/PAH source profiles_BaldwinEtAl2016_correctedUMO1.csv", stringsAsFactors=FALSE)

pmf2 <- merge(pmf, sources, by="Abbreviation", all.y=TRUE)

# create separate df for each pmf factor
f1profile <- select(pmf2, Abbreviation, f1_pctFactor,Power.plant.emissions,Residential.heating,Coal.average,
                    Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                    Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                    Asphalt,CT.T0, CT.T45, CT.T376,CT.dust.6,CT.dust.7)
f2profile <- select(pmf2, Abbreviation, f2_pctFactor,Power.plant.emissions,Residential.heating,Coal.average,
                    Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                    Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                    Asphalt,CT.T0, CT.T45, CT.T376,CT.dust.6,CT.dust.7)
f3profile <- select(pmf2, Abbreviation, f3_pctFactor,Power.plant.emissions,Residential.heating,Coal.average,
                    Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                    Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                    Asphalt,CT.T0, CT.T45, CT.T376,CT.dust.6,CT.dust.7)
f4profile <- select(pmf2, Abbreviation, f4_pctFactor,Power.plant.emissions,Residential.heating,Coal.average,
                    Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                    Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                    Asphalt,CT.T0, CT.T45, CT.T376,CT.dust.6,CT.dust.7)

# make each pmf factor df long
f1profileL <- melt(f1profile, id.vars=c("Abbreviation","f1_pctFactor"), variable.name = "source", value.name = "sourcePC")
# f1profileL$f1_pctFactor <- f1profileL$f1_pctFactor / 100

f2profileL <- melt(f2profile, id.vars=c("Abbreviation","f2_pctFactor"), variable.name = "source", value.name = "sourcePC")
# f2profileL$f2_pctFactor <- f2profileL$f2_pctFactor / 100

f3profileL <- melt(f3profile, id.vars=c("Abbreviation","f3_pctFactor"), variable.name = "source", value.name = "sourcePC")
# f3profileL$f3_pctFactor <- f3profileL$f3_pctFactor / 100

f4profileL <- melt(f4profile, id.vars=c("Abbreviation","f4_pctFactor"), variable.name = "source", value.name = "sourcePC")
# f4profileL$f4_pctFactor <- f4profileL$f4_pctFactor / 100


# plot all source profiles vs PMF factor 1
ggplot(f1profileL, aes(x=Abbreviation)) +
  geom_line(aes(y=sourcePC, group=source), color="red") +
  geom_point(aes(y=sourcePC),size=3, shape=16, color="red") +
  geom_line(aes(y=f1_pctFactor, group=1), color="black") +
  geom_point(aes(y=f1_pctFactor),size=3,shape=1, color="black") +
  facet_wrap(~ source, ncol=4) +
  xlab("Individual PAH compounds") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position="none", 
        axis.text.x = element_text(size=rel(.7), angle=90, hjust=1, vjust=.4),
        axis.text.y = element_text(size=rel(.7)),
        axis.title.y = element_text(size=rel(.7)),
        axis.title.x = element_text(size=rel(.7)))
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_glri11_Factor1.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_glri11_Factor1.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_glri11_Factor1.pdf"))

# plot all source profiles vs PMF factor 2
ggplot(f2profileL, aes(x=Abbreviation)) +
  geom_line(aes(y=sourcePC, group=source), color="red") +
  geom_point(aes(y=sourcePC),size=3, shape=16, color="red") +
  geom_line(aes(y=f2_pctFactor, group=1), color="black") +
  geom_point(aes(y=f2_pctFactor),size=3,shape=1, color="black") +
  facet_wrap(~ source, ncol=4) +
  xlab("Individual PAH compounds") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position="none", 
        axis.text.x = element_text(size=rel(.7), angle=90, hjust=1, vjust=.4),
        axis.text.y = element_text(size=rel(.7)),
        axis.title.y = element_text(size=rel(.7)),
        axis.title.x = element_text(size=rel(.7)))
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_glri11_Factor2.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_glri11_Factor2.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_glri11_Factor2.pdf"))

# plot all source profiles vs PMF factor 3
ggplot(f3profileL, aes(x=Abbreviation)) +
  geom_line(aes(y=sourcePC, group=source), color="red") +
  geom_point(aes(y=sourcePC),size=3, shape=16, color="red") +
  geom_line(aes(y=f3_pctFactor, group=1), color="black") +
  geom_point(aes(y=f3_pctFactor),size=3,shape=1, color="black") +
  facet_wrap(~ source, ncol=4) +
  xlab("Individual PAH compounds") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position="none", 
        axis.text.x = element_text(size=rel(.7), angle=90, hjust=1, vjust=.4),
        axis.text.y = element_text(size=rel(.7)),
        axis.title.y = element_text(size=rel(.7)),
        axis.title.x = element_text(size=rel(.7)))
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_glri10_Factor3.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_glri10_Factor3.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_glri10_Factor3.pdf"))

# plot all source profiles vs PMF factor 4
ggplot(f4profileL, aes(x=Abbreviation)) +
  geom_line(aes(y=sourcePC, group=source), color="red") +
  geom_point(aes(y=sourcePC),size=3, shape=16, color="red") +
  geom_line(aes(y=f4_pctFactor, group=1), color="black") +
  geom_point(aes(y=f4_pctFactor),size=3,shape=1, color="black") +
  facet_wrap(~ source, ncol=4) +
  xlab("Individual PAH compounds") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position="none", 
        axis.text.x = element_text(size=rel(.7), angle=90, hjust=1, vjust=.4),
        axis.text.y = element_text(size=rel(.7)),
        axis.title.y = element_text(size=rel(.7)),
        axis.title.x = element_text(size=rel(.7)))
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_TECsites_Factor4of4.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_TECsites_Factor4of4.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_TECsites_Factor4of4.pdf"))


#-----------------------------------------------------------------------
#  Compute chi-squared between proportional conc's 
#       PMF factor vs each source 
#   (PMF uncertainties not included - see uncertainties version below)
#-----------------------------------------------------------------------
# compute difference in mean of each PMF factor and source PC for each compound-
f1profileL$pcdiff <- abs(f1profileL$f1_pctFactor - f1profileL$sourcePC)
f2profileL$pcdiff <- abs(f2profileL$f2_pctFactor - f2profileL$sourcePC)
f3profileL$pcdiff <- abs(f3profileL$f3_pctFactor - f3profileL$sourcePC)
f4profileL$pcdiff <- abs(f4profileL$f4_pctFactor - f4profileL$sourcePC)

# chi2
f1profileL$pcdiff2 <- (f1profileL$pcdiff)^2
f1profileL$pcmean <- (f1profileL$f1_pctFactor + f1profileL$sourcePC)/2
f1profileL$chi2 <- f1profileL$pcdiff2 / f1profileL$pcmean
f1chi2summary <- group_by(f1profileL, source)%>%
  summarise(f1_chi2 = sum(chi2)) 

f2profileL$pcdiff2 <- (f2profileL$pcdiff)^2
f2profileL$pcmean <- (f2profileL$f2_pctFactor + f2profileL$sourcePC)/2
f2profileL$chi2 <- f2profileL$pcdiff2 / f2profileL$pcmean
f2chi2summary <- group_by(f2profileL, source)%>%
  summarise(f2_chi2 = sum(chi2))

f3profileL$pcdiff2 <- (f3profileL$pcdiff)^2
f3profileL$pcmean <- (f3profileL$f3_pctFactor + f3profileL$sourcePC)/2
f3profileL$chi2 <- f3profileL$pcdiff2 / f3profileL$pcmean
f3chi2summary <- group_by(f3profileL, source)%>%
  summarise(f3_chi2 = sum(chi2)) 

f4profileL$pcdiff2 <- (f4profileL$pcdiff)^2
f4profileL$pcmean <- (f4profileL$f4_pctFactor + f4profileL$sourcePC)/2
f4profileL$chi2 <- f4profileL$pcdiff2 / f4profileL$pcmean
f4chi2summary <- group_by(f4profileL, source)%>%
  summarise(f4_chi2 = sum(chi2)) 

# combine chi2 values for each pmf factor into one df
pmfChi2summary <- merge(f1chi2summary, f2chi2summary, by="source")
pmfChi2summary <- merge(pmfChi2summary, f3chi2summary, by="source")
pmfChi2summary <- merge(pmfChi2summary, f4chi2summary, by="source")

write.csv(pmfChi2summary, file = "C:/Users/akbaldwi/Documents/EPA PMF/Output/PMF_chi2Summary_glri11_2factor.csv", row.names = FALSE)


#-----------------------------------------------------------------------
#  Compute chi-squared between proportional conc's 
#       PMF factor vs each source 
#    ** PMF uncertainties included **
# * this is the method suggested by Gary Norris: Normalized Squared Difference
#-----------------------------------------------------------------------
# import csv with summary of PMF DISP uncertainties (DISP Avgs)
disp <- read.csv("C:/Users/akbaldwi/Documents/EPA PMF/Output/DISP_avgs_glri10and11.csv", stringsAsFactors=FALSE)
# delete all but the 12 compounds used in profiles
disp <- filter(disp, Abbrev != "")

# merge DISP values with profiles
f1profileL <- merge(f1profileL, disp, by.x="Abbreviation", by.y="Abbrev")
f2profileL <- merge(f2profileL, disp, by.x="Abbreviation", by.y="Abbrev")
f3profileL <- merge(f3profileL, disp, by.x="Abbreviation", by.y="Abbrev")
f4profileL <- merge(f4profileL, disp, by.x="Abbreviation", by.y="Abbrev")

# compute difference in mean of each PMF factor and source PC for each compound-
f1profileL$pcdiff <- abs(f1profileL$f1_pctFactor - f1profileL$sourcePC)
f2profileL$pcdiff <- abs(f2profileL$f2_pctFactor - f2profileL$sourcePC)
f3profileL$pcdiff <- abs(f3profileL$f3_pctFactor - f3profileL$sourcePC)
f4profileL$pcdiff <- abs(f4profileL$f4_pctFactor - f4profileL$sourcePC)

#----------------------------------------
# compute Normalized Squared Difference
#  ((sourcePC - samplePC) / DISP Error)^2

# glri11 - DISP Avg on Percent of Total
f1profileL$NSD_onPctTotal <- (f1profileL$pcdiff / f1profileL$DISPAvg_PctTotal_glri11_F1)^2
f2profileL$NSD_onPctTotal <- (f2profileL$pcdiff / f2profileL$DISPAvg_PctTotal_glri11_F2)^2
#-----
# glri10 - DISP Avg on Percent of Total
f1profileL$NSD_onPctTotal <- (f1profileL$pcdiff / f1profileL$DISPAvg_PctTotal_glri10_F1)^2
f2profileL$NSD_onPctTotal <- (f2profileL$pcdiff / f2profileL$DISPAvg_PctTotal_glri10_F2)^2
f3profileL$NSD_onPctTotal <- (f3profileL$pcdiff / f3profileL$DISPAvg_PctTotal_glri10_F3)^2
#-----
# Sum the Normalized Squared Differences for each source (equivalent to sum of chi2)
f1_sumNSD <- group_by(f1profileL, source)%>%
  summarize(f1_sumNSD = sum(NSD_onPctTotal))
f2_sumNSD <- group_by(f2profileL, source)%>%
  summarize(f2_sumNSD = sum(NSD_onPctTotal))
f3_sumNSD <- group_by(f3profileL, source)%>%
  summarize(f3_sumNSD = sum(NSD_onPctTotal))

# Merge sumNSD's of each factor into one dataframe
sumNSD <- merge(f1_sumNSD, f2_sumNSD, by="source")
sumNSD <- merge(sumNSD, f3_sumNSD, by="source")

write.csv(sumNSD, file = "C:/Users/akbaldwi/Documents/EPA PMF/Output/PMF_chi2_withUncertainty_Summary_glri11_2factor.csv", row.names = FALSE)
# OR
write.csv(sumNSD, file = "C:/Users/akbaldwi/Documents/EPA PMF/Output/PMF_chi2_withUncertainty_Summary_glri10_3factor.csv", row.names = FALSE)


#==============================================================
#       Concentrations of different PMF factors
#==============================================================
pmfConc <- select(pmf, Abbreviation, f1=f1_conc, f2=f2_conc)

pmfConc <- melt(pmfConc, id.vars="Abbreviation", variable.name = "factor", value.name = "concentration")

ggplot(pmfConc, aes(x=Abbreviation, y=concentration, color=factor))+
  geom_line(aes(group=factor))+
  ylab("Concentration, in ug/kg") +
  xlab("PAH compound") +
  theme_bw()+
  theme(legend.position=c(0.05,.95), legend.justification=c(0,1),
        legend.title=element_text(size=rel(1.)),
        legend.text=element_text(size=rel(1.)),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=rel(1.)),
        axis.text.y = element_text(size=rel(1.)),
        axis.title.y = element_text(size=rel(1.)),
        axis.title.x = element_text(size=rel(1.)))
ggsave(paste(pathToSave,"PMF_ConcMeans_glri11_factor.pdf"), width=4, height=3.5, dpi=300)
ggsave(paste(pathToSave,"PMF_ConcMeans_glri11_factor.png"), width=4, height=3.5, dpi=300)
shell.exec(paste(pathToSave,"PMF_ConcMeans_glri11_factor.pdf")) 












