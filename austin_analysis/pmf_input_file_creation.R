# Script to create the two input files needed for
# EPA's Positive Matrix Factorization (PMF) model.
# File 1 is sample concentrations in wide format,
# File 2 is uncertainties in wide format. 
# More info on these files in user guide EPA/600/R-14/108.
#------------------------------------------------------------

library(remake)
library(dplyr)
library(reshape2)


#=================================================================
# File 1: Sample concentrations
#=================================================================

#===========================================================
# OPTION 1: Create a PMF input file using just 12 compounds
#===========================================================
#   Use the dataframe created for the profiles analysis, 
#   but make it wide and merge compound names with casrn
profiles <- make('prepped_for_profiles')

# use the samples df to get casrn and param names
samples <- make('samples') # get sample data (no blanks or duplicates)
parms <- select(samples, param=PARAM_SYNONYM, casrn=CAS_NUM)
parms <- parms[!duplicated(parms),]  # remove duplicate rows

profiles <- merge(profiles, parms, by="casrn", all.x=TRUE) # merge profile data with param names

# Replace special characters in compound names
profiles$param <- gsub(",","", profiles$param)
profiles$param <- gsub("-","", profiles$param)
profiles$param <- gsub("(","_", profiles$param, fixed=TRUE)
profiles$param <- gsub(")","_", profiles$param)

# remove unwanted columns and make wide
profiles <- select(profiles, sample_id, param, RESULT)
profilesW <- dcast(profiles, sample_id ~ param)

# compute total conc per sample. Use to filter PMF data set
profilesW$Total <- rowSums(profilesW[,2:13])

# create input dataset of only samples exceeding TEC (1610 ppb)
profilesTEC <- filter(profilesW, Total > 1610) 

#===========================================================
# OPTION 2: Create a PMF input file using all compounds
#===========================================================
# nondetects:
#   concentration: 1/2 Detection Limit
#   uncertainty: 5/6 * DL  (as in Qi et al 2016, and PMF User Guide)

# note to Austin: need to include info on detect/nondetect, 
#   because don't want to include compounds with lots of nondetects

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

# create new result column, where result = 0.5 DL for nondetects
samplesAll$result2 <- ifelse(samplesAll$RESULT <= samplesAll$DETECT_LIMIT, 0.5*samplesAll$DETECT_LIMIT, samplesAll$RESULT)

# determine %nondetect for each compound
samplesAll$detect <- ifelse(samplesAll$RESULT <= samplesAll$DETECT_LIMIT, 0, 1)

nondetectSummary <- group_by(samplesAll, param) %>%
  summarise(n=n(), 
            totalDetects = sum(detect),
            pctDetect = (totalDetects/n)*100)

# omit compounds with <80% detects (based on summary above)
nondetectSummary <- filter(nondetectSummary, pctDetect > 80)
frequentCompounds <- nondetectSummary$param 

samplesAll <- filter(samplesAll, param %in% frequentCompounds)

# Replace special characters in compound names
samplesAll$param <- gsub(",","", samplesAll$param)
samplesAll$param <- gsub("-","", samplesAll$param)
samplesAll$param <- gsub("(","_", samplesAll$param, fixed=TRUE)
samplesAll$param <- gsub(")","_", samplesAll$param)

# remove unwanted columns 
samplesAll <- select(samplesAll, sample_id, param, result2)

# compute total sample conc
samplesAll <- group_by(samplesAll, sample_id) %>%
  mutate(totalConc = sum(result2))

# make wide
samplesAllW <- dcast(samplesAll, sample_id + totalConc ~ param, value.var = "result2")


#-----------------------------------------------------------
# File 2: Uncertainties
# For now use the stdev on duplicates as the uncertainty for each compound. 
#   Duplicate uncertainties are in qa_dulicates.
#-----------------------------------------------------------

#-------------------------------
# OPTION 1: 12-COMPOUNDS
#-------------------------------
qa_duplicates <- make('qa_duplicates')

samples12 <- filter(samples, sourceProfile12==TRUE)

# Merge samples with qa_duplicates
uncert <- merge(samples12, qa_duplicates, by="PARAM_SYNONYM")
uncert <- select(uncert, sample_id=unique_id, param=PARAM_SYNONYM, RESULT, DL= DETECT_LIMIT, mean, median, stdev)


# Compute uncertainty for each value using equation from Qi et al 16 (Nature; p.4)
#   uncert = sqrt((error x conc)^2 + DL^2)
#            for "error" use median RPD in duplicates 
uncert$uncertainty <- sqrt(((uncert$median / 100) * uncert$RESULT)^2 + uncert$DL^2)

uncert <- select(uncert, sample_id, param, uncertainty)

# Replace special characters in compound names
uncert$param <- gsub(",","", uncert$param)
uncert$param <- gsub("-","", uncert$param)

uncert$param <- gsub("(","_", uncert$param, fixed=TRUE)
uncert$param <- gsub(")","_", uncert$param)

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

uncertW <- merge(uncertW, total, by="sample_id")


#  uncertW and profilesW should have same number of rows, but don't. 
#  figure out which row to remove from uncertW (was likely omitted
#  from profilesW because it had a zero)
setdiff(uncertW$sample_id, profilesW$sample_id)

# remove sample MI-KAL from uncertW df - it's not in the profilesW df (must have a nondetect?)
uncertW <- filter(uncertW, sample_id != "MI-KAL")

# create uncertainty file with only samples exceeding TEC
tecSites <- profilesTEC$sample_id

uncertTEC <- filter(uncertW, sample_id %in% tecSites)

# write to csv for PMF input
write.csv(profilesW, file = "C:/Users/akbaldwi/Documents/EPA PMF/Data/concentrations_for_PMF.csv", row.names = FALSE)
write.csv(uncertW, file = "C:/Users/akbaldwi/Documents/EPA PMF/Data/uncertainties_for_PMF.csv", row.names = FALSE)

write.csv(profilesTEC, file = "C:/Users/akbaldwi/Documents/EPA PMF/Data/concentrationsTEC_for_PMF.csv", row.names = FALSE)
write.csv(uncertTEC, file = "C:/Users/akbaldwi/Documents/EPA PMF/Data/uncertaintiesTEC_for_PMF.csv", row.names = FALSE)

#----------------------------------------------------
#----------------------------------------------------
#-------------------------------
# OPTION 2: ALL COMPOUNDS with >80% detection
#-------------------------------
qa_duplicates <- make('qa_duplicates')

# Filter low detect compounds 
samples80 <- filter(samples, PARAM_SYNONYM %in% frequentCompounds)

# Merge samples with qa_duplicates
uncert <- merge(samples80, qa_duplicates, by="PARAM_SYNONYM")

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

# Replace special characters in compound names
uncert$param <- gsub(",","", uncert$param)
uncert$param <- gsub("-","", uncert$param)

uncert$param <- gsub("(","_", uncert$param, fixed=TRUE)
uncert$param <- gsub(")","_", uncert$param)

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

uncertW <- merge(uncertW, total, by="sample_id")


#  uncertW and profilesW should have same number of rows, but don't. 
#  figure out which row to remove from uncertW (was likely omitted
#  from profilesW because it had a zero)
setdiff(uncertW$sample_id, profilesW$sample_id)

# remove sample MI-KAL from uncertW df - it's not in the profilesW df (must have a nondetect?)
uncertW <- filter(uncertW, sample_id != "MI-KAL")

# create uncertainty file with only samples exceeding TEC
tecSites <- profilesTEC$sample_id

uncertTEC <- filter(uncertW, sample_id %in% tecSites)

# write to csv for PMF input
write.csv(profilesW, file = "C:/Users/akbaldwi/Documents/EPA PMF/Data/concentrations_for_PMF.csv", row.names = FALSE)
write.csv(uncertW, file = "C:/Users/akbaldwi/Documents/EPA PMF/Data/uncertainties_for_PMF.csv", row.names = FALSE)

write.csv(profilesTEC, file = "C:/Users/akbaldwi/Documents/EPA PMF/Data/concentrationsTEC_for_PMF.csv", row.names = FALSE)
write.csv(uncertTEC, file = "C:/Users/akbaldwi/Documents/EPA PMF/Data/uncertaintiesTEC_for_PMF.csv", row.names = FALSE)


