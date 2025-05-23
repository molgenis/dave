##################################################
# Dependencies (R equivalents of your Python libs)
##################################################
#install.packages(c("tidyverse","xgboost","caret","fastDummies","pROC"))
#install.packages("SHAPforxgboost") # if not installed

library(tidyverse)
library(xgboost)
library(caret)
library(fastDummies)
library(pROC)             # for AUC
library(crunch)     # To compress results
# library(SHAPforxgboost) # optional if you want further SHAP tools

#############################################
# Data import and split in LB/LP and VUS/CF #
#############################################

# Read CSV
freeze5 <- read.csv("/Users/joeri/git/vkgl-secretome-protein-stability/data/freeze5.csv.gz", stringsAsFactors = FALSE)

# Split into "labeled" set = LB or LP and "unlabeled" set = VUS or CF
freeze5_LB_LP <- freeze5 %>% filter(ann_classificationVKGL %in% c("LB", "LP"))
freeze5_VUS_CF <- freeze5 %>% filter(ann_classificationVKGL %in% c("VUS", "CF"))

##################################################
# Data preparation for ML on the LB/LP set
##################################################

# Select or remove features relevant for training
featureSelection <- c("ann_classificationVKGL", "delta_")
featureRemoval   <- c("delta_aaSeq", "delta_total.energy")

# Keep only columns whose names match any in featureSelection (a loose translation of regex) and drop columns from featureRemoval if they exist
freeze5_LB_LP_prep_for_ML <- freeze5_LB_LP %>%
  select(matches(paste(featureSelection, collapse="|"))) %>%
  select(-any_of(featureRemoval))

# Convert LB/LP labels into booleans
freeze5_LB_LP_prep_for_ML <- freeze5_LB_LP_prep_for_ML %>% mutate( ann_classificationVKGL = case_when(ann_classificationVKGL == "LB" ~ FALSE, ann_classificationVKGL == "LP" ~ TRUE), ann_classificationVKGL = as.logical(ann_classificationVKGL))

# Split into dependent (X) and independent (y)
X <- freeze5_LB_LP_prep_for_ML %>% select(-ann_classificationVKGL)
y <- freeze5_LB_LP_prep_for_ML$ann_classificationVKGL

# Random seed for reproducible splitting
set.seed(42)

# First split (80% train vs. 20% validation + test)
train_idx <- createDataPartition(y, p = 0.8, list = FALSE)
X_train   <- X[train_idx, ]
y_train   <- y[train_idx]
X_valtest <- X[-train_idx, ]
y_valtest <- y[-train_idx]

# Second split (50% validation vs. 50% test)
val_idx <- createDataPartition(y_valtest, p = 0.5, list = FALSE)
X_val   <- X_valtest[val_idx, ]
y_val   <- y_valtest[val_idx]
X_test  <- X_valtest[-val_idx, ]
y_test  <- y_valtest[-val_idx]

##############################
# Machine learning procedure #
##############################

# Create xgb.DMatrix objects
dtrain <- xgb.DMatrix(as.matrix(X_train), label = as.numeric(y_train))
dval   <- xgb.DMatrix(as.matrix(X_val),   label = as.numeric(y_val))
dtest  <- xgb.DMatrix(as.matrix(X_test),  label = as.numeric(y_test))

# Define hyperparameters
params <- list(objective = "binary:logistic", eval_metric = "logloss", eta = 0.01, subsample = 0.6, learning_rate = 0.09, seed = 42)

# Train with watchlist (train vs. val)
watchlist <- list(train = dtrain, validation = dval)
xgb_reg <- xgb.train(params = params, data = dtrain, nrounds = 150, watchlist = watchlist, verbose = 1)

# Check AUC on validation
val_pred <- predict(xgb_reg, newdata = dval)
val_auc  <- auc(y_val, val_pred)
cat("XGBoost model validation AUC:", val_auc, "\n")

# Re-train final model by combining train + val
dtrain_val <- xgb.DMatrix(as.matrix(rbind(X_train, X_val)), label = as.numeric(c(y_train, y_val)))
xgb_reg_final <- xgb.train(params = params, data = dtrain_val, nrounds = 150, verbose = 1)

# Evaluate on test set (which was held out until now)
test_pred <- predict(xgb_reg_final, newdata = as.matrix(X_test))
test_auc  <- auc(y_test, test_pred)
cat("XGBoost model final test AUC:", test_auc, "\n")

##################################################
# Explain predictions (SHAP approach)
##################################################

# Define logistic functions
sigmoid <- function(x) { 1 / (1 + exp(-x)) }
inverseSigmoid <- function(x) { log(x / (1 - x)) }

# Function that emulates your Python "predict_as_prob_scaled_SHAP"
#    Using XGBoost’s predcontrib=TRUE to get SHAP (log-odds) contributions + bias term.
predict_as_prob_scaled_SHAP <- function(df, model, split_name) {
  # Turn df into matrix
  df_mat <- as.matrix(df)
  
  # SHAP values (in log-odds); last column is the BIAS (base value)
  shap_contrib <- predict(model, newdata = df_mat, predcontrib = TRUE)
  
  # Separate out the base_value (bias) from the per-feature SHAP
  base_value <- shap_contrib[, ncol(shap_contrib)]
  shap_values <- shap_contrib[, -ncol(shap_contrib), drop = FALSE]
  
  # Convert base_value to probabilities
  base_prob <- sigmoid(base_value)
  
  # Sum of all SHAP log-odds for each row
  shap_sum <- rowSums(shap_values)
  final_log_odds <- base_value + shap_sum
  final_prob <- sigmoid(final_log_odds)
  
  # Scale the per-feature SHAP by fraction of total shap_sum
  # (Add a small epsilon if shap_sum is zero to avoid dividing by zero.)
  eps <- 1e-15
  shap_values_scaled <- shap_values / (shap_sum + eps)
  
  # Probability difference from base to final
  prob_diff <- final_prob - base_prob
  
  # Multiply scaled shap by the probability difference
  shap_values_scaled_prob <- shap_values_scaled * prob_diff
  
  # Build output dataframe
  # Rename feature columns with ".sph"
  colnames(shap_values_scaled_prob) <- paste0(colnames(df), ".sph")
  result <- as.data.frame(shap_values_scaled_prob)
  
  # Add columns for base and final probabilities
  result$BaseProbability.sph  <- base_prob
  result$FinalProbability.sph <- final_prob
  result$MLSplit              <- split_name
  
  return(result)
}

# 3) Apply the function to get SHAP-based scaled probabilities
X_train_SHAP_prob_values <- predict_as_prob_scaled_SHAP(X_train, xgb_reg_final, "Train")
X_val_SHAP_prob_values   <- predict_as_prob_scaled_SHAP(X_val,   xgb_reg_final, "Train") 
X_test_SHAP_prob_values  <- predict_as_prob_scaled_SHAP(X_test,  xgb_reg_final, "Test")

##################################################
# Prepare the VUS+CF set, get predictions
##################################################
fr5_VUS_CF_forPred <- freeze5_VUS_CF %>%
  # Keep the columns you used for training
  select(matches(paste(featureSelection, collapse = "|"))) %>%
  # Drop the same featureRemoval columns
  select(-any_of(featureRemoval)) %>%
  # Drop 'ann_classificationVKGL' if present
  select(-any_of("ann_classificationVKGL"))

fr5_VUS_CF_SHAP_prob_values <- predict_as_prob_scaled_SHAP(
  fr5_VUS_CF_forPred, xgb_reg_final, "Unknown"
)

##################################################
# Combine all results with original data, save
##################################################
allPredsInOne <- bind_rows(
  X_train_SHAP_prob_values,
  X_val_SHAP_prob_values,
  X_test_SHAP_prob_values,
  fr5_VUS_CF_SHAP_prob_values
)

freeze5_plus_pred <- bind_cols(freeze5, allPredsInOne)

write.csv.gz(freeze5_plus_pred, "/Users/joeri/git/vkgl-secretome-protein-stability/data/freeze5_predictions_R.csv.gz", row.names = FALSE, quote = FALSE)


