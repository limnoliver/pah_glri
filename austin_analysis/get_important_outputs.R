# Documentation of various steps in the PAH analysis
# and how/where to get the appropriate data and figures
# you may need.
library(remake)
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


# 3. Remake everything by running
make()

# 4. 'COMMIT' AND 'PUSH' YOUR WORK
# when working in pah_glri (say, adding your own scripts in 'austin_analysis'), 
# 'commit' your changes (with a commit message describing what you did) periodically 
# (I do this about as often as a hit the save button -- be sure to hit save, then commit)
# 'push' changes (using the up arrow in the git tab) periodically but less frequently than 
# committing. I usually do this if I walk away from my computer, or if I've done a major change/
# addition. This won't make changes to my repository (limnoliver/pah_glri) but rather
# to your own fork (akbaldwi-usgs/pah_glri).



##################################################
######### get sample concentrations ##############
# raw data are retrieved (1_get_raw_data), processed (2_process_data), 
# and filtered (3_filter_data) to get to "samples"

# target "samples" is all of the concentration data
# that has been merged with compound-specific data
# (e.g., low vs high molecular weight). BDLs have been
# converted to zero, and raw concentration data has been merged
# with compound-specific data (e.g., molecular weight) using function ?pah::get_compound_info
# compound metadata stored in pah package pah::pah_compounds

samples <- make('samples') # just samples, no blanks or duplicates
# each row is a different sample-compound combination (so a long table)
# units are ppb
# site and compound metadata columns at end of table
pah16 <- make('prepped_totals') # sample ID's and corresponding EPA16 concentration
site_order <- make('sample_order_pah16') # a list of site ID's, from lowest to highest EPA16 concentration

##################################################
######### QA data ##################

# these targets were made in 3_filter_data
blanks <- make('blanks') # blanks with coresponding samples and duplicates
duplicates <- make('duplicates') # duplicates with corresponding samples

# summary QA targets were made in 10_qa_summary

# for each compound, blank summary statistics were calculated:
# n_FB_bdl: the number of instances where the blank was below detection limit
# n_both_bdl: the number of times the blank and corresponding sample were BDL
# n_SA_bdl: the number of times the corresponding sample was BDL, but not the blank
# n_both_adl: the number of times both the blank and corresponding sample were above
# detection limit. These were the instances that were analyzed for percent differences
# mean: mean of all percent difference between blank and sample (when both were > DL)
# medain: as above
# stdev: as above
# min: as above
# max: as above
qa_blanks <- make('qa_blanks') 

# for each duplicate sample, summary statistics were calculated relative to the sample.
# relative percent differences were calculated for each sample-compound combination, 
# and then for each compound, summary statistics were calculated:
# mean: for each compound, mean of all RPD between duplicate and sample
# median, stdev, min, max: as above
# n_both_bdl: number of samples where both duplicate and sample were BDL
# n_one_bdl: number of samples where one but not both sample and duplicate were BDL
qa_duplicates <- make('qa_duplicates')

##################################################
######### profiles ###############################

# step 6_profile_data

# profile dat, list of 2
# part 1: long format for profiles, where each row is a unique site-compound-source combination
# and both sample and source proportional concentration have their own columns
# part 2: sum chi2 for each sample-source combination

profiles <- make('profiles')

# this is the step where you can also include or drop creosote as a source, which
# only has info for 11 compounds. The above target is without creosote. To get the 
# profiles with creosote:

profiles_creosote <- make('profiles_creosote')

##################################################
######### PCA anaysis ############################

# step 8_pca_analysis

# summary stats of PCA are output when the pah::pah_pca
# function is run, but can also be accessed in the output
# pah_pca lets you change selection criteria for PCA components,
# but in this case, we kept all components that explained >= 10% variation

pca_dat <- make('pca_no_creosote')

# PCA stats output
pca_dat$pca_summary # first 4 components all explained > 10%, together explained 78% of total variance

# raw PCA output of selected components
# output in long format, where each row is a sample or source, and column "type" 
# tells you whether it is a sample or source value
raw_pca <- pca_dat$pca_dat

# pca distances, long format (each row is unique sample-source combo)
# gives euclidean distances in all selected PCA component space
pca_dist <- pca_dat$pca_distance

# closest source by site using PCA distance
top_pca <- make('pca_top_sources_bysite')

##################################################
######### Double ratios ##########################

# step 7_double_ratios

# raw ratios for source and samples
# uses pah::source_ratios and pah::calc_ratios
# Each ratio is its own column, and sources and samples are both in "sample_id"
# column but are distinguished by the "sample_type" column.
ratios <- make('ratios')

# calculated distances for each double ratio plot using pah::cal_ratio_dist
# list of 3: raw (raw distances for each source-sample-double ratio plot combo), source (raw distance
# data summarized by source), and samples (top sources for each sample calculated in each double ratio plot)
# see 7_double_ratios.yml for more descriptions
ratio_dist <- make('ratio_distance')

##################################################
######### parent vs alkyl and ESBTU ##############

# step 9_parent_weight

# the compound metadata to calculate these are already in 
# the "samples" target and come from pah::pah_compounds
# The data could be retrieved by using pah::pah_mw_parent function and 
# setting plot = FALSE which will output a data frame
# see 9_parent_weight.yml for how this function is used

##################################################
######### Percent mass fractions ####

# step 11_mass_fractions
# get percent mass fraction based on summary statistics of all samples
mf_summary <- make('percent_by_weight_summary')

# get percent mass fraction for all sample-source combinations
mf_bysample <- make('percent_by_weight_bysample')



