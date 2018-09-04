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
pmf <- read.csv("C:/Users/akbaldwi/Documents/EPA PMF/Output/glri_2vs3_defVsMod_summary_for_R.csv", stringsAsFactors=FALSE)

# import source profiles from literature
sources <- read.csv("C:/Users/akbaldwi/Documents/MMSD_PAHs/R/for_GLRI_reuse/PAH source profiles_BaldwinEtAl2016_correctedUMO1.csv", stringsAsFactors=FALSE)

pmf2 <- merge(pmf, sources, by="Abbreviation", all.y=TRUE)

# create separate df for each pmf factor
f1_glri2def <- select(pmf2, Abbreviation, f1_glri2def=f1_pctFactor_glri2def,Power.plant.emissions,Residential.heating,Coal.average,
                    Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                    Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                    Asphalt,CT.T0, CT.T45, CT.T376,CT.dust.6)
f2_glri2def <- select(pmf2, Abbreviation, f2_glri2def=f2_pctFactor_glri2def,Power.plant.emissions,Residential.heating,Coal.average,
                      Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                      Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                      Asphalt,CT.T0, CT.T45, CT.T376,CT.dust.6)
f1_glri3def <- select(pmf2, Abbreviation, f1_glri3def=f1_pctFactor_glri3def,Power.plant.emissions,Residential.heating,Coal.average,
                      Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                      Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                      Asphalt,CT.T0, CT.T45, CT.T376,CT.dust.6)
f2_glri3def <- select(pmf2, Abbreviation, f2_glri3def=f2_pctFactor_glri3def,Power.plant.emissions,Residential.heating,Coal.average,
                      Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                      Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                      Asphalt,CT.T0, CT.T45, CT.T376,CT.dust.6)
f3_glri3def <- select(pmf2, Abbreviation, f3_glri3def=f3_pctFactor_glri3def,Power.plant.emissions,Residential.heating,Coal.average,
                      Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                      Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                      Asphalt,CT.T0, CT.T45, CT.T376,CT.dust.6)
f1_glri2mod <- select(pmf2, Abbreviation, f1_glri2mod=f1_pctFactor_glri2mod,Power.plant.emissions,Residential.heating,Coal.average,
                      Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                      Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                      Asphalt,CT.T0, CT.T45, CT.T376,CT.dust.6)
f2_glri2mod <- select(pmf2, Abbreviation, f2_glri2mod=f2_pctFactor_glri2mod,Power.plant.emissions,Residential.heating,Coal.average,
                      Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                      Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                      Asphalt,CT.T0, CT.T45, CT.T376,CT.dust.6)
f1_glri3mod <- select(pmf2, Abbreviation, f1_glri3mod=f1_pctFactor_glri3mod,Power.plant.emissions,Residential.heating,Coal.average,
                      Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                      Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                      Asphalt,CT.T0, CT.T45, CT.T376,CT.dust.6)
f2_glri3mod <- select(pmf2, Abbreviation, f2_glri3mod=f2_pctFactor_glri3mod,Power.plant.emissions,Residential.heating,Coal.average,
                      Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                      Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                      Asphalt,CT.T0, CT.T45, CT.T376,CT.dust.6)
f3_glri3mod <- select(pmf2, Abbreviation, f3_glri3mod=f3_pctFactor_glri3mod,Power.plant.emissions,Residential.heating,Coal.average,
                      Coke.oven.emissions,Diesel.vehicle,Gasoline.vehicle,Traffic.tunnel.air,Vehicle.traffic.avg,Used.motor.oil.1,
                      Used.motor.oil.2,Pine.combustion.1,Pine.combustion.2,Oak.combustion,Fuel.oil.combustion,Tire.particles,
                      Asphalt,CT.T0, CT.T45, CT.T376,CT.dust.6)


# make each pmf factor df long
f1_glri2defL <- melt(f1_glri2def, id.vars=c("Abbreviation","f1_glri2def"), variable.name = "source", value.name = "sourcePC")
f2_glri2defL <- melt(f2_glri2def, id.vars=c("Abbreviation","f2_glri2def"), variable.name = "source", value.name = "sourcePC")
f1_glri3defL <- melt(f1_glri3def, id.vars=c("Abbreviation","f1_glri3def"), variable.name = "source", value.name = "sourcePC")
f2_glri3defL <- melt(f2_glri3def, id.vars=c("Abbreviation","f2_glri3def"), variable.name = "source", value.name = "sourcePC")
f3_glri3defL <- melt(f3_glri3def, id.vars=c("Abbreviation","f3_glri3def"), variable.name = "source", value.name = "sourcePC")
f1_glri2modL <- melt(f1_glri2mod, id.vars=c("Abbreviation","f1_glri2mod"), variable.name = "source", value.name = "sourcePC")
f2_glri2modL <- melt(f2_glri2mod, id.vars=c("Abbreviation","f2_glri2mod"), variable.name = "source", value.name = "sourcePC")
f1_glri3modL <- melt(f1_glri3mod, id.vars=c("Abbreviation","f1_glri3mod"), variable.name = "source", value.name = "sourcePC")
f2_glri3modL <- melt(f2_glri3mod, id.vars=c("Abbreviation","f2_glri3mod"), variable.name = "source", value.name = "sourcePC")
f3_glri3modL <- melt(f3_glri3mod, id.vars=c("Abbreviation","f3_glri3mod"), variable.name = "source", value.name = "sourcePC")


# plot all source profiles vs each PMF factor
ggplot(f1_glri2defL, aes(x=Abbreviation)) +
  geom_line(aes(y=sourcePC, group=source), color="red") +
  geom_point(aes(y=sourcePC),size=3, shape=16, color="red") +
  geom_line(aes(y=f1_glri2def, group=1), color="black") +
  geom_point(aes(y=f1_glri2def),size=3,shape=1, color="black") +
  facet_wrap(~ source, ncol=4) +
  xlab("Individual PAH compounds") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position="none", 
        axis.text.x = element_text(size=rel(.7), angle=90, hjust=1, vjust=.4),
        axis.text.y = element_text(size=rel(.7)),
        axis.title.y = element_text(size=rel(.7)),
        axis.title.x = element_text(size=rel(.7)))
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f1_glri2def.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f1_glri2def.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f1_glri2def.pdf"))

ggplot(f2_glri2defL, aes(x=Abbreviation)) +
  geom_line(aes(y=sourcePC, group=source), color="red") +
  geom_point(aes(y=sourcePC),size=3, shape=16, color="red") +
  geom_line(aes(y=f2_glri2def, group=1), color="black") +
  geom_point(aes(y=f2_glri2def),size=3,shape=1, color="black") +
  facet_wrap(~ source, ncol=4) +
  xlab("Individual PAH compounds") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position="none", 
        axis.text.x = element_text(size=rel(.7), angle=90, hjust=1, vjust=.4),
        axis.text.y = element_text(size=rel(.7)),
        axis.title.y = element_text(size=rel(.7)),
        axis.title.x = element_text(size=rel(.7)))
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f2_glri2def.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f2_glri2def.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f2_glri2def.pdf"))

ggplot(f1_glri3defL, aes(x=Abbreviation)) +
  geom_line(aes(y=sourcePC, group=source), color="red") +
  geom_point(aes(y=sourcePC),size=3, shape=16, color="red") +
  geom_line(aes(y=f1_glri3def, group=1), color="black") +
  geom_point(aes(y=f1_glri3def),size=3,shape=1, color="black") +
  facet_wrap(~ source, ncol=4) +
  xlab("Individual PAH compounds") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position="none", 
        axis.text.x = element_text(size=rel(.7), angle=90, hjust=1, vjust=.4),
        axis.text.y = element_text(size=rel(.7)),
        axis.title.y = element_text(size=rel(.7)),
        axis.title.x = element_text(size=rel(.7)))
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f1_glri3def.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f1_glri3def.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f1_glri3def.pdf"))

ggplot(f2_glri3defL, aes(x=Abbreviation)) +
  geom_line(aes(y=sourcePC, group=source), color="red") +
  geom_point(aes(y=sourcePC),size=3, shape=16, color="red") +
  geom_line(aes(y=f2_glri3def, group=1), color="black") +
  geom_point(aes(y=f2_glri3def),size=3,shape=1, color="black") +
  facet_wrap(~ source, ncol=4) +
  xlab("Individual PAH compounds") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position="none", 
        axis.text.x = element_text(size=rel(.7), angle=90, hjust=1, vjust=.4),
        axis.text.y = element_text(size=rel(.7)),
        axis.title.y = element_text(size=rel(.7)),
        axis.title.x = element_text(size=rel(.7)))
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f2_glri3def.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f2_glri3def.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f2_glri3def.pdf"))

ggplot(f3_glri3defL, aes(x=Abbreviation)) +
  geom_line(aes(y=sourcePC, group=source), color="red") +
  geom_point(aes(y=sourcePC),size=3, shape=16, color="red") +
  geom_line(aes(y=f3_glri3def, group=1), color="black") +
  geom_point(aes(y=f3_glri3def),size=3,shape=1, color="black") +
  facet_wrap(~ source, ncol=4) +
  xlab("Individual PAH compounds") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position="none", 
        axis.text.x = element_text(size=rel(.7), angle=90, hjust=1, vjust=.4),
        axis.text.y = element_text(size=rel(.7)),
        axis.title.y = element_text(size=rel(.7)),
        axis.title.x = element_text(size=rel(.7)))
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f3_glri3def.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f3_glri3def.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f3_glri3def.pdf"))

ggplot(f1_glri2modL, aes(x=Abbreviation)) +
  geom_line(aes(y=sourcePC, group=source), color="red") +
  geom_point(aes(y=sourcePC),size=3, shape=16, color="red") +
  geom_line(aes(y=f1_glri2mod, group=1), color="black") +
  geom_point(aes(y=f1_glri2mod),size=3,shape=1, color="black") +
  facet_wrap(~ source, ncol=4) +
  xlab("Individual PAH compounds") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position="none", 
        axis.text.x = element_text(size=rel(.7), angle=90, hjust=1, vjust=.4),
        axis.text.y = element_text(size=rel(.7)),
        axis.title.y = element_text(size=rel(.7)),
        axis.title.x = element_text(size=rel(.7)))
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f1_glri2mod.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f1_glri2mod.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f1_glri2mod.pdf"))

ggplot(f2_glri2modL, aes(x=Abbreviation)) +
  geom_line(aes(y=sourcePC, group=source), color="red") +
  geom_point(aes(y=sourcePC),size=3, shape=16, color="red") +
  geom_line(aes(y=f2_glri2mod, group=1), color="black") +
  geom_point(aes(y=f2_glri2mod),size=3,shape=1, color="black") +
  facet_wrap(~ source, ncol=4) +
  xlab("Individual PAH compounds") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position="none", 
        axis.text.x = element_text(size=rel(.7), angle=90, hjust=1, vjust=.4),
        axis.text.y = element_text(size=rel(.7)),
        axis.title.y = element_text(size=rel(.7)),
        axis.title.x = element_text(size=rel(.7)))
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f2_glri2mod.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f2_glri2mod.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f2_glri2mod.pdf"))

ggplot(f1_glri3modL, aes(x=Abbreviation)) +
  geom_line(aes(y=sourcePC, group=source), color="red") +
  geom_point(aes(y=sourcePC),size=3, shape=16, color="red") +
  geom_line(aes(y=f1_glri3mod, group=1), color="black") +
  geom_point(aes(y=f1_glri3mod),size=3,shape=1, color="black") +
  facet_wrap(~ source, ncol=4) +
  xlab("Individual PAH compounds") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position="none", 
        axis.text.x = element_text(size=rel(.7), angle=90, hjust=1, vjust=.4),
        axis.text.y = element_text(size=rel(.7)),
        axis.title.y = element_text(size=rel(.7)),
        axis.title.x = element_text(size=rel(.7)))
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f1_glri3mod.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f1_glri3mod.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f1_glri3mod.pdf"))

ggplot(f2_glri3modL, aes(x=Abbreviation)) +
  geom_line(aes(y=sourcePC, group=source), color="red") +
  geom_point(aes(y=sourcePC),size=3, shape=16, color="red") +
  geom_line(aes(y=f2_glri3mod, group=1), color="black") +
  geom_point(aes(y=f2_glri3mod),size=3,shape=1, color="black") +
  facet_wrap(~ source, ncol=4) +
  xlab("Individual PAH compounds") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position="none", 
        axis.text.x = element_text(size=rel(.7), angle=90, hjust=1, vjust=.4),
        axis.text.y = element_text(size=rel(.7)),
        axis.title.y = element_text(size=rel(.7)),
        axis.title.x = element_text(size=rel(.7)))
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f2_glri3mod.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f2_glri3mod.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f2_glri3mod.pdf"))

ggplot(f3_glri3modL, aes(x=Abbreviation)) +
  geom_line(aes(y=sourcePC, group=source), color="red") +
  geom_point(aes(y=sourcePC),size=3, shape=16, color="red") +
  geom_line(aes(y=f3_glri3mod, group=1), color="black") +
  geom_point(aes(y=f3_glri3mod),size=3,shape=1, color="black") +
  facet_wrap(~ source, ncol=4) +
  xlab("Individual PAH compounds") +
  ylab("PAH proportional concentrations") +
  theme_bw()  + 
  theme(legend.position="none", 
        axis.text.x = element_text(size=rel(.7), angle=90, hjust=1, vjust=.4),
        axis.text.y = element_text(size=rel(.7)),
        axis.title.y = element_text(size=rel(.7)),
        axis.title.x = element_text(size=rel(.7)))
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f3_glri3mod.pdf"), width=11, height=8, dpi=300)
ggsave(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f3_glri3mod.png"), width=11, height=8, dpi=300)
shell.exec(paste(pathToSave,"PAH_profiles_AllSources_vs_PMF_f3_glri3mod.pdf"))



#-----------------------------------------------------------------------
#  Compute chi-squared between proportional conc's 
#       PMF factor vs each source 
#   (PMF uncertainties not included - see uncertainties version below)
#-----------------------------------------------------------------------
# compute difference in mean of each PMF factor and source PC for each compound-
f1_glri2defL$pcdiff <- abs(f1_glri2defL$f1_glri2def - f1_glri2defL$sourcePC)
f2_glri2defL$pcdiff <- abs(f2_glri2defL$f2_glri2def - f2_glri2defL$sourcePC)
f1_glri3defL$pcdiff <- abs(f1_glri3defL$f1_glri3def - f1_glri3defL$sourcePC)
f2_glri3defL$pcdiff <- abs(f2_glri3defL$f2_glri3def - f2_glri3defL$sourcePC)
f3_glri3defL$pcdiff <- abs(f3_glri3defL$f3_glri3def - f3_glri3defL$sourcePC)
f1_glri2modL$pcdiff <- abs(f1_glri2modL$f1_glri2mod - f1_glri2modL$sourcePC)
f2_glri2modL$pcdiff <- abs(f2_glri2modL$f2_glri2mod - f2_glri2modL$sourcePC)
f1_glri3modL$pcdiff <- abs(f1_glri3modL$f1_glri3mod - f1_glri3modL$sourcePC)
f2_glri3modL$pcdiff <- abs(f2_glri3modL$f2_glri3mod - f2_glri3modL$sourcePC)
f3_glri3modL$pcdiff <- abs(f3_glri3modL$f3_glri3mod - f3_glri3modL$sourcePC)


# chi2
f1_glri2defL$pcdiff2 <- (f1_glri2defL$pcdiff)^2
f1_glri2defL$pcmean <- (f1_glri2defL$f1_glri2def + f1_glri2defL$sourcePC)/2
f1_glri2defL$chi2 <- f1_glri2defL$pcdiff2 / f1_glri2defL$pcmean
f1_glri2def_chi2summary <- group_by(f1_glri2defL, source)%>%
  summarise(f1_glri2def = sum(chi2)) 

f2_glri2defL$pcdiff2 <- (f2_glri2defL$pcdiff)^2
f2_glri2defL$pcmean <- (f2_glri2defL$f2_glri2def + f2_glri2defL$sourcePC)/2
f2_glri2defL$chi2 <- f2_glri2defL$pcdiff2 / f2_glri2defL$pcmean
f2_glri2def_chi2summary <- group_by(f2_glri2defL, source)%>%
  summarise(f2_glri2def = sum(chi2)) 

f1_glri3defL$pcdiff2 <- (f1_glri3defL$pcdiff)^2
f1_glri3defL$pcmean <- (f1_glri3defL$f1_glri3def + f1_glri3defL$sourcePC)/2
f1_glri3defL$chi2 <- f1_glri3defL$pcdiff2 / f1_glri3defL$pcmean
f1_glri3def_chi2summary <- group_by(f1_glri3defL, source)%>%
  summarise(f1_glri3def = sum(chi2)) 

f2_glri3defL$pcdiff2 <- (f2_glri3defL$pcdiff)^2
f2_glri3defL$pcmean <- (f2_glri3defL$f2_glri3def + f2_glri3defL$sourcePC)/2
f2_glri3defL$chi2 <- f2_glri3defL$pcdiff2 / f2_glri3defL$pcmean
f2_glri3def_chi2summary <- group_by(f2_glri3defL, source)%>%
  summarise(f2_glri3def = sum(chi2)) 

f3_glri3defL$pcdiff2 <- (f3_glri3defL$pcdiff)^2
f3_glri3defL$pcmean <- (f3_glri3defL$f3_glri3def + f3_glri3defL$sourcePC)/2
f3_glri3defL$chi2 <- f3_glri3defL$pcdiff2 / f3_glri3defL$pcmean
f3_glri3def_chi2summary <- group_by(f3_glri3defL, source)%>%
  summarise(f3_glri3def = sum(chi2)) 

f1_glri2modL$pcdiff2 <- (f1_glri2modL$pcdiff)^2
f1_glri2modL$pcmean <- (f1_glri2modL$f1_glri2mod + f1_glri2modL$sourcePC)/2
f1_glri2modL$chi2 <- f1_glri2modL$pcdiff2 / f1_glri2modL$pcmean
f1_glri2mod_chi2summary <- group_by(f1_glri2modL, source)%>%
  summarise(f1_glri2mod = sum(chi2)) 

f2_glri2modL$pcdiff2 <- (f2_glri2modL$pcdiff)^2
f2_glri2modL$pcmean <- (f2_glri2modL$f2_glri2mod + f2_glri2modL$sourcePC)/2
f2_glri2modL$chi2 <- f2_glri2modL$pcdiff2 / f2_glri2modL$pcmean
f2_glri2mod_chi2summary <- group_by(f2_glri2modL, source)%>%
  summarise(f2_glri2mod = sum(chi2)) 

f1_glri3modL$pcdiff2 <- (f1_glri3modL$pcdiff)^2
f1_glri3modL$pcmean <- (f1_glri3modL$f1_glri3mod + f1_glri3modL$sourcePC)/2
f1_glri3modL$chi2 <- f1_glri3modL$pcdiff2 / f1_glri3modL$pcmean
f1_glri3mod_chi2summary <- group_by(f1_glri3modL, source)%>%
  summarise(f1_glri3mod = sum(chi2)) 

f2_glri3modL$pcdiff2 <- (f2_glri3modL$pcdiff)^2
f2_glri3modL$pcmean <- (f2_glri3modL$f2_glri3mod + f2_glri3modL$sourcePC)/2
f2_glri3modL$chi2 <- f2_glri3modL$pcdiff2 / f2_glri3modL$pcmean
f2_glri3mod_chi2summary <- group_by(f2_glri3modL, source)%>%
  summarise(f2_glri3mod = sum(chi2)) 

f3_glri3modL$pcdiff2 <- (f3_glri3modL$pcdiff)^2
f3_glri3modL$pcmean <- (f3_glri3modL$f3_glri3mod + f3_glri3modL$sourcePC)/2
f3_glri3modL$chi2 <- f3_glri3modL$pcdiff2 / f3_glri3modL$pcmean
f3_glri3mod_chi2summary <- group_by(f3_glri3modL, source)%>%
  summarise(f3_glri3mod = sum(chi2)) 


# combine chi2 values for each pmf factor into one df
pmfChi2summary <- merge(f1_glri2def_chi2summary, f2_glri2def_chi2summary, by="source")
pmfChi2summary <- merge(pmfChi2summary, f1_glri3def_chi2summary, by="source")
pmfChi2summary <- merge(pmfChi2summary, f2_glri3def_chi2summary, by="source")
pmfChi2summary <- merge(pmfChi2summary, f3_glri3def_chi2summary, by="source")
pmfChi2summary <- merge(pmfChi2summary, f1_glri2mod_chi2summary, by="source")
pmfChi2summary <- merge(pmfChi2summary, f2_glri2mod_chi2summary, by="source")
pmfChi2summary <- merge(pmfChi2summary, f1_glri3mod_chi2summary, by="source")
pmfChi2summary <- merge(pmfChi2summary, f2_glri3mod_chi2summary, by="source")
pmfChi2summary <- merge(pmfChi2summary, f3_glri3mod_chi2summary, by="source")

write.csv(pmfChi2summary, file = "C:/Users/akbaldwi/Documents/EPA PMF/Output/PMF_chi2Summary_2v3factor_and_convergenceCriteriaComparison.csv", row.names = FALSE)


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












