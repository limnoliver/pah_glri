# Script to compare PMF factors to one another, to determine whether they're
# different or basically the same. The 3-factor solution has two factors that are fairly
# similar; this script looks at their similarity, and how they compare to the 
# corresponding factor in the 2-factor solution. 
# This script also explores whether the solutions are different using the default vs modified
# convergence criteria. 
#=========================================================================================

pathToSave <- paste("C:/Users/akbaldwi/Documents/EPA PMF/Output/", sep="")

library(reshape2)
library(ggplot2)
library(dplyr)
library(matrixStats)
library(RColorBrewer)
library(Hmisc)
#====================================================================

# import PMF factors for 2-factor solution and 3-factor solution
#   and using default vs modified convergence criteria
pmf <- read.csv("C:/Users/akbaldwi/Documents/EPA PMF/Output/glri_2vs3_defVsMod_summary_for_R.csv", stringsAsFactors=FALSE)

#==========================================================================
# Compute Spearman correlations: rcorr method (Hmisc package)
#==========================================================================

# Correlation matrix: 
corr <- select(pmf, 2:11)
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
fileToSave <- paste(pathToSave, "Correlation_matrix_PMF_2vs3Factors_and_DefaultVsModifiedConvergenceCriteria.csv",sep="/")
write.table(profileCorrSpearman_Corr, fileToSave, row.names=FALSE, sep=",")
fileToSave <- paste(pathToSave, "Correlation_matrix_PMF_2vs3Factors_and_DefaultVsModifiedConvergenceCriteria_PVALUES.csv",sep="/")
write.table(profileCorrSpearman_Pval, fileToSave, row.names=FALSE, sep=",")

