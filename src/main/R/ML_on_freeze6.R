library(crunch)     # To compress results
library(R.utils)    # for 'gunzip', 'mkdirs'
library(dplyr)      # Data manipulation
library(tidyr)      # Data manipulation
library(xgboost)
library(caret)
library(randomForest)


#######################
# Adjustable settings #
#######################
rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability" # root directory that contains README.md, data/, img/, out/, src/, etc.
ftThr <- 5 # minimum nr of both LB/LP variants when considering a sequence feature

###################################
# Invariants: seperators and seed #
###################################
CONCAT_SEP <- "~" # used to combine UniProt hierarchy: category -> type -> description -> etc
VALUE_SEP <- "|" # used to combine multiple annotation values in the output data freeze
set.seed(42)

##################################
# Load data and helper functions #
##################################
frz6loc <- paste(rootDir, "data", "freeze6.csv.gz", sep="/")
frz6 <- read.csv(frz6loc)
source(paste(rootDir, "src", "main", "R", "ML_on_freeze6_functions.R", sep="/"))


####
# Prep, 
####
# Cleanup of variables, dummify etc
prepFrz6 <- prepForML(frz6)
# Before we do anything else, keep 20% of data unused in a balanced way
train_idx <- createDataPartition(prepFrz6$ann_classificationVKGL, p = 0.8, list = FALSE)
prepFrz6_train <- prepFrz6[train_idx,] # use for training
prepFrz6_test <- prepFrz6[-train_idx,] # do not touch until the very end
# Find out which sequence features are usable
usableSeqFt <- seqFtWithMinimumNrPerLabel(prepFrz6_train, ftThr)
cat(paste0(length(usableSeqFt), " sequence features have ", ftThr, " (or more) benign and ", ftThr, " (or more) pathogenic variants" ))
# Iterate over usable sequence features and train micromodels
for(sf in usableSeqFt){
  sf <- "DOMAINS_AND_SITES" # debug
  cat(paste0("Creating micromodel for ", sf))
  
  sfDat <- sliceBySeqFt(prepFrz6_train, sf)
  model <- trainModelOn(sfDat)
  
}




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

