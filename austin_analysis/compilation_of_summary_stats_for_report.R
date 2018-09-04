# Script to compile summary stats from different PAH source ID analyses
#  for discussion in report

library(remake)
library(tidyr)
library(reshape2)
devtools::install_github('limnoliver/pah')
samples <- make('samples') # just samples, no blanks or duplicates
pah16 <- make('prepped_totals') # sample ID's and corresponding EPA16 concentration
pah16 <- select(pah16, -Priority16_bin)

#=================================
#        Parent / Alkyl
#=================================
pa <- drop_na(samples, parentAlkyl)

# what compounds are included in parents? alkyls?
parents <- filter(pa, parentAlkyl == "parent")
unique(parents$PARAM_SYNONYM)
alkyls <- filter(pa, parentAlkyl == "alkyl")
unique(alkyls$PARAM_SYNONYM)
unique(samples$PARAM_SYNONYM)

paSummary <- group_by(pa, unique_id, parentAlkyl) %>%
  summarise(n=n(),
            sumConc = sum(RESULT, na.rm=TRUE))

paSummary <- select(paSummary, -n)

paSummary <- dcast(paSummary, unique_id ~ parentAlkyl, value.var = "sumConc")

paSummary$parentAlkylRatio <- paSummary$parent / paSummary$alkyl

paSummary <- select(paSummary, unique_id, parentAlkylRatio)

#=================================
#        HMW / LMW
#=================================
mw <- drop_na(samples, molwt_highlow)

mwSummary <- group_by(mw, unique_id, molwt_highlow) %>%
  summarise(n=n(),
            sumConc = sum(RESULT, na.rm=TRUE))

mwSummary <- select(mwSummary, -n)

mwSummary <- dcast(mwSummary, unique_id ~ molwt_highlow, value.var = "sumConc")

mwSummary$molwtRatio <- mwSummary$HMW / mwSummary$LMW

mwSummary <- select(mwSummary, unique_id, molwtRatio)


#=================================
#        TOC
#=================================
toc <- filter(samples, PARAM_SYNONYM == "TOC")
toc <- select(toc, unique_id, TOC=RESULT)
 
#=================================
# Combine above into single table
#=================================
sampleSummary <- merge(toc, pah16, by.x="unique_id", by.y="sample_id", all=TRUE)
sampleSummary <- merge(sampleSummary, paSummary, by="unique_id", all=TRUE)
sampleSummary <- merge(sampleSummary, mwSummary, by="unique_id", all=TRUE)

# Next... merge PECQ, TECQ, and sumESBTU into sampleSummary 


write.csv(sampleSummary, file = "C:/Users/akbaldwi/Documents/GLRI/PAHs/Manuscript/sample_concs_etc_summary_table.csv", row.names = FALSE)
