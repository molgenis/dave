library(crunch)     # To compress results
library(R.utils)    # for 'gunzip', 'mkdirs'
library(dplyr)      # Data manipulation
library(tidyr)      # Data manipulation
library(xgboost)
library(caret)
library(randomForest)
library(reshape2)
library(pheatmap)
library(rstanarm)

#######################
# Adjustable settings #
#######################
rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability" # root directory that contains README.md, data/, img/, out/, src/, etc.
ftThr <- 10 # minimum nr of both LB/LP variants when considering a sequence feature
options(mc.cores = parallel::detectCores())


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


##########################################################
# Clean, slice, investigate features and train submodels #
##########################################################
# Cleanup of variables, dummify etc
prepFrz6 <- prepForML(frz6)
# Before we do anything else, keep 20% of data unused in a balanced way
train_idx <- createDataPartition(prepFrz6$ann_classificationVKGL, p = 0.8, list = FALSE)
prepFrz6_train <- prepFrz6[train_idx,] # use for training
prepFrz6_test <- prepFrz6[-train_idx,] # do not touch until the very end
# Find out which sequence features are usable
usableSeqFt <- seqFtWithMinimumNrPerLabel(prepFrz6_train, ftThr)
cat(paste0(length(usableSeqFt), " sequence features have ", ftThr, " (or more) benign and ", ftThr, " (or more) pathogenic variants" ))
# Iterate over usable sequence features and train submodels
all_est <- data.frame()
for(sf in usableSeqFt){
  #sf <- "DOMAINS_AND_SITES~DOMAIN~Ig-like C2-type" # debug
  cat(paste0("Running BayesGLM for ", sf, "\n"))
  sfDat <- sliceBySeqFt(prepFrz6_train, sf)
  model <- trainModelOn(sfDat)
  coefs <- coef(model)
  sf_coef_est <- data.frame(
    sf = sf,
    coefficient = names(coefs),
    estimate = as.numeric(coefs),
    row.names = NULL
  )
  all_est <- rbind(all_est, sf_coef_est)
}
# Transform to log odds ratios to probabilities via 1 / (1 + exp(-estimate))
all_est <- all_est %>% mutate(prob_scale = plogis(estimate))
# Convert into wide format and replace NA with 0
coef_wide <- dcast(all_est, sf ~ coefficient, value.var = "estimate") # estimate or prob_scale
coef_wide[is.na(coef_wide)] <- 0
# Show as heatmap for inspection
coef_wide$`(Intercept)` <- NULL # leave out Intercept ?
mat <- as.matrix(coef_wide[,-1])
rownames(mat) <- coef_wide$sf
pheatmap(mat, color = rev(viridis(100)), cluster_rows = TRUE, cluster_cols = TRUE, scale = "none") # plasma, magma


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

