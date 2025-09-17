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


############################################
# Define which features will be target, 'explainers', or removed.
#######################################
# TODO: better engineering, e.g. cumu_bin is similar to cumu_prob, remove more?
# prediction target
topLvl <- "ann_classificationVKGL"
# "explainers" that will be linked to feature importance
midLvl <- c("delta_total.energy",
            "delta_ligand_nr_of_predicted_pockets",
            "delta_DNAb_binding_affinity_pKd",
            "delta_RNAb_binding_affinity_pKd",
            "delta_DNAs_cumu_bin",
            "delta_RNAs_cumu_bin",
            "delta_ProtS_cumu_bin",
            "abs_delta_total.energy",
            "abs_delta_ligand_nr_of_predicted_pockets",
            "abs_delta_DNAb_binding_affinity_pKd",
            "abs_delta_RNAb_binding_affinity_pKd",
            "abs_delta_DNAs_cumu_bin",
            "abs_delta_RNAs_cumu_bin",
            "abs_delta_ProtS_cumu_bin"
            )
# not usable for training, remove these features
removeFeat <- c("gene",
                "TranscriptID",
                "UniProtID",
                "dna_variant_chrom",
                "dna_variant_pos",
                "dna_variant_ref",
                "dna_variant_alt",
                "dna_variant_assembly",
                "ann_proteinIschaperoned",
                "ann_proteinLocalization",
                "delta_aaSeq",
                "seqFt",
                "chain")


##########################################################
# Clean, slice, investigate features and train submodels #
##########################################################
# Cleanup of variables, dummify etc
prepFrz6 <- commonPrepForAllPred(frz6, removeFeat)
# Before we do anything else, keep 20% of data unused in a balanced way
train_idx <- createDataPartition(prepFrz6$ann_classificationVKGL, p = 0.8, list = FALSE)
prepFrz6_train <- prepFrz6[train_idx,] # use for training
prepFrz6_test <- prepFrz6[-train_idx,] # do not touch until the very end


# iterate over prediction targets and make models

for(target in midLvl){
  #target <- "delta_total.energy" # debug
  cat(paste0("Training model for ", target, "\n"))
  # remove all midLvl and topLvl, but not the current prediction target
  prepFrz6_train2 <- prepFrz6_train %>% select(-all_of(c(setdiff(midLvl, target), topLvl)))
  
  x <- prepFrz6_train2 %>% select(-all_of(target))
  y <- prepFrz6_train2 %>% pull(all_of(target))
  rf_model <- randomForest(x = x, y = y)
  
  preds <- predict(rf_model, newdata = prepFrz6_test)
  
  # R²
  rsq <- 1 - sum(((prepFrz6_test[[target]]) - preds)^2) / sum((prepFrz6_test[[target]] - mean(prepFrz6_test[[target]]))^2)
  rsq
  cat(paste0("rsq is ", rsq, "\n"))
  
  
  
  #rf_model <- randomForest(target ~ ., data = prepFrz6_train2)
  #model <- trainModel(prepFrz6_train, target)
}

# variable importance
importance(rf_model)

# nice plot
varImpPlot(rf_model)


topLvlTrain <- prepFrz6_train %>% select(-any_of(c(midLvl, removeFeat)))

#topLvlTrain <- prepFrz6_train %>% select(matches(paste(featureSelection, collapse="|"))) %>% select(-any_of(c(midLvl, removeThis)))


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

