# Script to create the two input files needed for
# EPA's Positive Matrix Factorization (PMF) model.
# File 1 is sample concentrations in wide format,
# File 2 is uncertainties in wide format. 
# More info on these files in user guide EPA/600/R-14/108.

# Include all compounds with detection frequencies > 80%
#------------------------------------------------------------

library(remake)
library(dplyr)
library(reshape2)

qa_duplicates <- make('qa_duplicates') 

#=================================================================
# File 1: Sample concentrations
#=================================================================

# nondetects:
#   concentration: 1/2 Detection Limit
#   uncertainty: 5/6 * DL  (as in Qi et al 2016, and PMF User Guide)

#   First determine %nondetect for each compound
#   Omit compounds with less than 80% detection

samples <- make('samples') # get sample data (no blanks or duplicates)
samplesAll <- select(samples, sample_id=unique_id, param=PARAM_SYNONYM, RESULT, DETECT_LIMIT)

# create list of params that aren't PAHs and omit them
nonPAH <- c("TOTAL_SOLIDS","TOC","GRAVEL_MED_4.75MM","GRAVEL_FINE_2MM","SAND_VCOARSE_0.85MM",
            "SAND_COARSE_0.425MM","SAND_MED_0.25MM","SAND_FINE_0.106MM","SAND_VFINE_0.075MM",
            "75.0_UM","31.3_UM","15.6_UM","7.8_UM","3.9_UM","1.95_UM","0.98_UM","Total PAH",
            "Priority Pollutant PAH","Acenaphthene-d10","Benzo(a)pyrene-d12","Naphthalene-d8",
            "Phenanthrene-d10","PCT_MOIST")
samplesAll <- filter(samplesAll, !param %in% nonPAH)

# also filter to only those params included in the "duplicates" df
params_qa <- qa_duplicates$PARAM_SYNONYM
samplesAll <- filter(samplesAll, param %in% params_qa)

# create new result column, where result = 0.5 DL for nondetects
samplesAll$result2 <- ifelse(samplesAll$RESULT <= samplesAll$DETECT_LIMIT, 0.5*samplesAll$DETECT_LIMIT, samplesAll$RESULT)

# determine %nondetect for each COMPOUND
samplesAll$detect <- ifelse(samplesAll$RESULT <= samplesAll$DETECT_LIMIT, 0, 1)

nondetectSummary <- group_by(samplesAll, param) %>%
  summarise(n=n(), 
            totalDetects = sum(detect),
            pctDetect = (totalDetects/n)*100)

# write detection summary to csv
write.csv(nondetectSummary, file = "C:/Users/akbaldwi/Documents/EPA PMF/Data/compound_detection_frequencies_by_compound.csv", row.names = FALSE)

# determine %nondetect for each SAMPLE
nondetectSummaryBySample <- group_by(samplesAll, sample_id) %>%
  summarise(n=n(), 
            totalDetects = sum(detect),
            pctDetect = (totalDetects/n)*100)

# write detection summary to csv
write.csv(nondetectSummaryBySample, file = "C:/Users/akbaldwi/Documents/EPA PMF/Data/compound_detection_frequencies_by_sample.csv", row.names = FALSE)


#-------- Not running this section on advice from Pentti -------------
# omit compounds with <80% detects (based on summary above)
nondetectSummary <- filter(nondetectSummary, pctDetect > 80)
frequentCompounds <- nondetectSummary$param 
samplesAll <- filter(samplesAll, param %in% frequentCompounds)
#---------------------------------------------------------------------

# Replace special characters in compound names
samplesAll$param <- gsub(",","", samplesAll$param)
samplesAll$param <- gsub("-","", samplesAll$param)
samplesAll$param <- gsub("(","_", samplesAll$param, fixed=TRUE)
samplesAll$param <- gsub(")","_", samplesAll$param)

# remove unwanted columns 
samplesAll <- select(samplesAll, sample_id, param, result2)

# compute total sample conc
samplesAll <- group_by(samplesAll, sample_id) %>%
  mutate(Total = sum(result2))

# make wide
samplesAllW <- dcast(samplesAll, sample_id + Total ~ param, value.var = "result2")


#-----------------------------------------------------------
# File 2: Uncertainties
# For now use the stdev on duplicates as the uncertainty for each compound. 
#   Duplicate uncertainties are in qa_dulicates.
#-----------------------------------------------------------

qa_duplicates <- make('qa_duplicates')

samplesU <- samples

# Replace special characters in compound names
samplesU$PARAM_SYNONYM <- gsub(",","", samplesU$PARAM_SYNONYM)
samplesU$PARAM_SYNONYM <- gsub("-","", samplesU$PARAM_SYNONYM)
samplesU$PARAM_SYNONYM <- gsub("(","_", samplesU$PARAM_SYNONYM, fixed=TRUE)
samplesU$PARAM_SYNONYM <- gsub(")","_", samplesU$PARAM_SYNONYM)

qa_duplicates$PARAM_SYNONYM <- gsub(",","", qa_duplicates$PARAM_SYNONYM)
qa_duplicates$PARAM_SYNONYM <- gsub("-","", qa_duplicates$PARAM_SYNONYM)
qa_duplicates$PARAM_SYNONYM <- gsub("(","_", qa_duplicates$PARAM_SYNONYM, fixed=TRUE)
qa_duplicates$PARAM_SYNONYM <- gsub(")","_", qa_duplicates$PARAM_SYNONYM)

# Filter to compounds in samplesAll df above
keepCompounds <- unique(samplesAll$param)
samplesU <- filter(samplesU, PARAM_SYNONYM %in% keepCompounds)

# Merge samples with qa_duplicates
uncert <- merge(samplesU, qa_duplicates, by="PARAM_SYNONYM", all.x=TRUE)

# create new result column, where result = 0.5 DL for nondetects
uncert$result2 <- ifelse(uncert$RESULT <= uncert$DETECT_LIMIT, 0.5*uncert$DETECT_LIMIT, uncert$RESULT)

uncert <- select(uncert, sample_id=unique_id, param=PARAM_SYNONYM, result2, DL= DETECT_LIMIT, mean, median, stdev)

# Compute uncertainty for each value using equation from Qi et al 16 (Nature; p.4)
#   Detections:
#         uncert = sqrt((error x conc)^2 + DL^2)
#            for "error" use median RPD in duplicates 
#   Nondetections:
#         uncert = 5/6 DL
uncert$uncertDetect <- sqrt(((uncert$median / 100) * uncert$result2)^2 + uncert$DL^2)
uncert$uncertNondetect <- uncert$DL * (5/6)

uncert$uncertainty <- ifelse(uncert$result2 < uncert$DL, uncert$uncertNondetect, uncert$uncertDetect)

uncert <- select(uncert, sample_id, param, uncertainty)

#======================================
# compute combined (total) uncertainty
#======================================
# combine uncertainties of indiv compounds using root sum of the squares
uncert$squared <- uncert$uncertainty^2

# sum squares 
uncert <- group_by(uncert, sample_id) %>%
  mutate(sumOfSquares = sum(squared))

# sqrt of sum of squares
uncert$rootSumOfSquares <- sqrt(uncert$sumOfSquares)

total <- select(uncert, sample_id, Total=rootSumOfSquares)

total <- unique(total)
#======================================

# make wide and merge Total column
uncertW <- select(uncert, sample_id, param, uncertainty)

uncertW <- dcast(uncertW, sample_id ~ param)

uncertW <- merge(total, uncertW, by="sample_id")

#=========================================
# write to csv for PMF input
write.csv(samplesAllW, file = "C:/Users/akbaldwi/Documents/EPA PMF/Data/concentrations_for_PMF_allCompounds.csv", row.names = FALSE)
write.csv(uncertW, file = "C:/Users/akbaldwi/Documents/EPA PMF/Data/uncertainties_for_PMF_allCompounds.csv", row.names = FALSE)

#=========================================
# create and export summary of nondetects per SAMPLE





