#  Script to help determine whether the uncertainties in the PMF input file 
#    are appropriate, or whether they need to be adjusted. 
#    (Based on email guidance from Norris and Pentti)

library(remake)
library(dplyr)
library(reshape2)
library(ggplot2)

conc <- read.csv("C:/Users/akbaldwi/Documents/EPA PMF/Data/concentrations_for_PMF.csv", stringsAsFactors=FALSE)
resids <- read.csv("C:/Users/akbaldwi/Documents/EPA PMF/Output/glri03_residuals_forR.csv", stringsAsFactors=FALSE)
scaledResids <- read.csv("C:/Users/akbaldwi/Documents/EPA PMF/Output/glri03_scaledresiduals_forR.csv", stringsAsFactors=FALSE)

#======================================================================================

# make dataframes long and merge
concW <- melt(conc, id.vars="sample_id", variable.name = "compound", value.name = "concentration")
concW <- rename(concW, Sample_IDs=sample_id)
residsW <- melt(resids, id.vars=c("Sample_IDs","Base_Run"), variable.name = "compound", value.name = "residuals")
scaledResidsW <- melt(scaledResids, id.vars=c("Sample_IDs","Base_Run"), variable.name = "compound", value.name = "scaledResiduals")

conres <- merge(concW, residsW, by=c("Sample_IDs","compound"))
conres <- merge(conres, scaledResidsW, by=c("Sample_IDs","compound","Base_Run"))

# compute fitted value (concentration + residual)
conres$fittedValue <- conres$concentration + conres$residuals

# plot fitted values vs scaled residuals
ggplot(conres, aes(x=fittedValue, y=scaledResiduals)) +
  geom_point() +
  facet_wrap(~ compound) +
  scale_x_log10()+
  labs(x="Fitted values\n(measured + residual)", y="Scaled residuals")+
  theme_bw()+
  theme(strip.text.x = element_text(size=rel(1)),
        strip.text.y = element_text(size=rel(1)),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10),
        axis.title.x = element_text(size=rel(1.)),
        axis.title.y = element_text(size=rel(1.)))
ggsave("C:/Users/akbaldwi/Documents/EPA PMF/Data/PMF_conc_uncertainty_assessment.pdf", width=8, height=8, dpi=300)
ggsave(("C:/Users/akbaldwi/Documents/EPA PMF/Data/PMF_conc_uncertainty_assessment.png"), width=8, height=8, dpi=300)
shell.exec("C:/Users/akbaldwi/Documents/EPA PMF/Data/PMF_conc_uncertainty_assessment.pdf")
                