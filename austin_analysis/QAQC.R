# Austin's script for QA-related summaries and plots. 
#  - detection limits
#  - blank results
#  - replicate results
########################################################

library(remake)
library(dplyr)
library(modeest)
##################################################
######### GENERAL WORFLOW REMINDERS ########

# 1. GET LATEST CHANGES TO pah_glri
# before you start working, it's a good idea to pull the latest version of both the 
# PAH package and this repository. From the Git tab > Shell, 'git pull upstream master' 
# will pull my most recent changes to the "pah_glri" repository. 

# 2. GET LATEST VERSION OF pah PACKAGE
# For now, since you don't have R tools installed, you can just pull the PAH package 
# and use it like you would any other package. Update it by running from the 
# RStudio console:
devtools::install_github('limnoliver/pah')

samples <- make('samples') # just samples, no blanks or duplicates

#=====================================================
#   detection limits
#=====================================================
# by compound
dl <- filter(samples, CLASS == "PAH") %>%
  group_by(DESCR) %>%
  summarize(min = min(DETECT_LIMIT),
            median = median(DETECT_LIMIT),
            mode  = mlv(DETECT_LIMIT, method='mfv')[['M']],
            max = max(DETECT_LIMIT))
  
# across all compounds
dl2 <- filter(samples, CLASS == "PAH") 
summary(dl2$DETECT_LIMIT)
  
  
  
  
  
  
  
  
  
  