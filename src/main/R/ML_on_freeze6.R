library(crunch)     # To compress results
library(R.utils)    # for 'gunzip', 'mkdirs'
library(dplyr)      # Data manipulation
library(tidyr)      # Data manipulation

#######################
# Adjustable settings #
#######################
rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability" # root directory that contains README.md, data/, img/, out/, src/, etc.


########################
# Seperator characters #
########################
CONCAT_SEP <- "~" # used to combine UniProt hierarchy: category -> type -> description -> etc
VALUE_SEP <- "|" # used to combine multiple annotation values in the output data freeze


#############
# Load data #
#############
frz6loc <- paste(rootDir, "data", "freeze6.csv.gz", sep="/")
frz6 <- read.csv(frz6loc)

##
# Train many 'micro models' on subsets of data using protein sequence features
##


##
# Train complete model on standard 80/20 split, all features and/or micromodels
##


####
# Apply complete model on VUS variants to get predictions
####


##
# Explain predictions using SHAP heuristic
##


#####
# Combine all results and save to final predictions and explanations
####

