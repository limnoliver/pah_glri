# Script to compare PMF output profiles to source profiles using plots and chi-squared.
#=========================================================================================

pathToSave <- paste("C:/Users/akbaldwi/Documents/GLRI/PAHs/PMF/Pentti_files/", sep="")

library(reshape2)
library(ggplot2)
library(dplyr)
library(matrixStats)
library(RColorBrewer)
#====================================================================

# import PMF factors
pmf <- read.csv("C:/Users/akbaldwi/Documents/GLRI/PAHs/PMF/Pentti_files/sedimF3FF3_source_profile_factors.csv", stringsAsFactors=FALSE)
# import source profiles from literature
sources <- read.csv("C:/Users/akbaldwi/Documents/MMSD_PAHs/R/for_GLRI_reuse/PAH source profiles_BaldwinEtAl2016_correctedUMO1.csv", stringsAsFactors=FALSE)
# import CTDust weathering profiles from Van Metre and Mahler volatilization study
ctWeathering <- read.csv("C:/Users/akbaldwi/Documents/MMSD_PAHs/Response to OReilly/weathering.AKB.proportionalConcsOnly.csv")
ctWeathering <- select(ctWeathering, -PAH, -CTD6, -CTD7)

pmf2 <- merge(pmf, sources, by="Abbreviation", all.y=TRUE)
pmf2 <- merge(pmf2, ctWeathering, by.x="Abbreviation", by.y="abbrev", all=TRUE)

# create separate df for each pmf factor
f1profile <- select(pmf2, Abbreviation, f1_pctFactor,Power.plant.emissions,Residential.heating,Coal.average,
                    Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                    Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                    Asphalt,CT.dust.6,CT.dust.7,CT_T0,CT_T1,CT_T5,CT_T45,CT_T149,CT_T232,CT_T328,CT_T376)
f2profile <- select(pmf2, Abbreviation, f2_pctFactor,Power.plant.emissions,Residential.heating,Coal.average,
                    Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                    Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                    Asphalt,CT.dust.6,CT.dust.7,CT_T0,CT_T1,CT_T5,CT_T45,CT_T149,CT_T232,CT_T328,CT_T376)
f3profile <- select(pmf2, Abbreviation, f3_pctFactor,Power.plant.emissions,Residential.heating,Coal.average,
                    Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                    Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                    Asphalt,CT.dust.6,CT.dust.7,CT_T0,CT_T1,CT_T5,CT_T45,CT_T149,CT_T232,CT_T328,CT_T376)
f4profile <- select(pmf2, Abbreviation, f4_pctFactor,Power.plant.emissions,Residential.heating,Coal.average,
                    Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                    Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                    Asphalt,CT.dust.6,CT.dust.7,CT_T0,CT_T1,CT_T5,CT_T45,CT_T149,CT_T232,CT_T328,CT_T376)
f5profile <- select(pmf2, Abbreviation, f5_pctFactor,Power.plant.emissions,Residential.heating,Coal.average,
                    Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                    Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                    Asphalt,CT.dust.6,CT.dust.7,CT_T0,CT_T1,CT_T5,CT_T45,CT_T149,CT_T232,CT_T328,CT_T376)
f6profile <- select(pmf2, Abbreviation, f6_pctFactor,Power.plant.emissions,Residential.heating,Coal.average,
                    Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                    Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                    Asphalt,CT.dust.6,CT.dust.7,CT_T0,CT_T1,CT_T5,CT_T45,CT_T149,CT_T232,CT_T328,CT_T376)

# make each pmf factor df long
f1profileL <- melt(f1profile, id.vars=c("Abbreviation","f1_pctFactor"), variable.name = "source", value.name = "sourcePC")
# f1profileL$f1_pctFactor <- f1profileL$f1_pctFactor / 100

f2profileL <- melt(f2profile, id.vars=c("Abbreviation","f2_pctFactor"), variable.name = "source", value.name = "sourcePC")
# f2profileL$f2_pctFactor <- f2profileL$f2_pctFactor / 100

f3profileL <- melt(f3profile, id.vars=c("Abbreviation","f3_pctFactor"), variable.name = "source", value.name = "sourcePC")
# f3profileL$f3_pctFactor <- f3profileL$f3_pctFactor / 100

f4profileL <- melt(f4profile, id.vars=c("Abbreviation","f4_pctFactor"), variable.name = "source", value.name = "sourcePC")
# f4profileL$f4_pctFactor <- f4profileL$f4_pctFactor / 100

f5profileL <- melt(f5profile, id.vars=c("Abbreviation","f5_pctFactor"), variable.name = "source", value.name = "sourcePC")

f6profileL <- melt(f6profile, id.vars=c("Abbreviation","f6_pctFactor"), variable.name = "source", value.name = "sourcePC")


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
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor1.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor1.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor1.pdf"))

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
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor2.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor2.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor2.pdf"))

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
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor3.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor3.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor3.pdf"))

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
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor4.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor4.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor4.pdf"))

# plot all source profiles vs PMF factor 5
ggplot(f5profileL, aes(x=Abbreviation)) +
  geom_line(aes(y=sourcePC, group=source), color="red") +
  geom_point(aes(y=sourcePC),size=3, shape=16, color="red") +
  geom_line(aes(y=f5_pctFactor, group=1), color="black") +
  geom_point(aes(y=f5_pctFactor),size=3,shape=1, color="black") +
  facet_wrap(~ source, ncol=4) +
  xlab("Individual PAH compounds") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position="none", 
        axis.text.x = element_text(size=rel(.7), angle=90, hjust=1, vjust=.4),
        axis.text.y = element_text(size=rel(.7)),
        axis.title.y = element_text(size=rel(.7)),
        axis.title.x = element_text(size=rel(.7)))
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor5.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor5.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor5.pdf"))

# plot all source profiles vs PMF factor 6
ggplot(f6profileL, aes(x=Abbreviation)) +
  geom_line(aes(y=sourcePC, group=source), color="red") +
  geom_point(aes(y=sourcePC),size=3, shape=16, color="red") +
  geom_line(aes(y=f6_pctFactor, group=1), color="black") +
  geom_point(aes(y=f6_pctFactor),size=3,shape=1, color="black") +
  facet_wrap(~ source, ncol=4) +
  xlab("Individual PAH compounds") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position="none", 
        axis.text.x = element_text(size=rel(.7), angle=90, hjust=1, vjust=.4),
        axis.text.y = element_text(size=rel(.7)),
        axis.title.y = element_text(size=rel(.7)),
        axis.title.x = element_text(size=rel(.7)))
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor6.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor6.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_Pentti_6factor_Factor6.pdf"))


#-----------------------------------------------------------------------
#  Compute chi-squared between proportional conc's 
#       PMF factor vs each source 
#-----------------------------------------------------------------------
# compute difference in mean of each PMF factor and source PC for each compound-
f1profileL$pcdiff <- abs(f1profileL$f1_pctFactor - f1profileL$sourcePC)
f2profileL$pcdiff <- abs(f2profileL$f2_pctFactor - f2profileL$sourcePC)
f3profileL$pcdiff <- abs(f3profileL$f3_pctFactor - f3profileL$sourcePC)
f4profileL$pcdiff <- abs(f4profileL$f4_pctFactor - f4profileL$sourcePC)
f5profileL$pcdiff <- abs(f5profileL$f5_pctFactor - f5profileL$sourcePC)
f6profileL$pcdiff <- abs(f6profileL$f6_pctFactor - f6profileL$sourcePC)

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

f5profileL$pcdiff2 <- (f5profileL$pcdiff)^2
f5profileL$pcmean <- (f5profileL$f5_pctFactor + f5profileL$sourcePC)/2
f5profileL$chi2 <- f5profileL$pcdiff2 / f5profileL$pcmean
f5chi2summary <- group_by(f5profileL, source)%>%
  summarise(f5_chi2 = sum(chi2)) 

f6profileL$pcdiff2 <- (f6profileL$pcdiff)^2
f6profileL$pcmean <- (f6profileL$f6_pctFactor + f6profileL$sourcePC)/2
f6profileL$chi2 <- f6profileL$pcdiff2 / f6profileL$pcmean
f6chi2summary <- group_by(f6profileL, source)%>%
  summarise(f6_chi2 = sum(chi2)) 

# combine chi2 values for each pmf factor into one df
pmfChi2summary <- merge(f1chi2summary, f2chi2summary, by="source")
pmfChi2summary <- merge(pmfChi2summary, f3chi2summary, by="source")
pmfChi2summary <- merge(pmfChi2summary, f4chi2summary, by="source")
pmfChi2summary <- merge(pmfChi2summary, f5chi2summary, by="source")
pmfChi2summary <- merge(pmfChi2summary, f6chi2summary, by="source")

write.csv(pmfChi2summary, file = "C:/Users/akbaldwi/Documents/GLRI/PAHs/PMF/Pentti_files/PMF_chi2Summary_Pentti_F3FF3.csv", row.names = FALSE)



#==============================================================
#       Concentrations of different PMF factors
#==============================================================
pmfConc <- select(pmf, Abbreviation, f1=f1_conc, f2=f2_conc, f3=f3_conc)

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
ggsave(paste(pathToSave,"PMF_ConcMeans_glri10_factor.pdf"), width=4, height=3.5, dpi=300)
ggsave(paste(pathToSave,"PMF_ConcMeans_glri10_factor.png"), width=4, height=3.5, dpi=300)
shell.exec(paste(pathToSave,"PMF_ConcMeans_glri10_factor.pdf")) 












