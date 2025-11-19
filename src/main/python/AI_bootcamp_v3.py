# Dependencies
###############
# Using Python 3.12.4
import pandas as pd
import xgboost as xgb
import numpy as np
import shap
import math
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score


# Data import and initial split
###############################
# Import data
freeze5 = pd.read_csv("/Users/joeri/git/dave/data/freeze5-provisional-NAs.csv.gz", low_memory=False)
# Split into on "labeled" test/validation/train data and "unlabeled" data that we want to apply the model on
freeze5_LB_LP = freeze5[(freeze5['ann_classificationVKGL'] == 'LB') | (freeze5['ann_classificationVKGL'] == 'LP')]
freeze5_VUS_CF = freeze5[(freeze5['ann_classificationVKGL'] == 'VUS') | (freeze5['ann_classificationVKGL'] == 'CF')]


# Data preparation for ML on the LB/LP set
##########################################
# Change protein localization from 3 categoricals into 2 booleans: "xx_membrane" and
# "xx_secreted" with the 3rd one ("xx_intracellular") being implicit if the other ones are false
freeze5_LB_LP_prep_for_ML = pd.get_dummies(freeze5_LB_LP, columns=['ann_proteinLocalization'], drop_first=True)
# Select columns that are relevant features for training the model
featureSelection = ["ann_classificationVKGL", "delta_"] #  "ann_proteinLocalization", "ann_proteinIschaperoned" --> not adding these because they're not functional!
featureRemoval = ["delta_aaSeq", "delta_total.energy"]
freeze5_LB_LP_prep_for_ML = freeze5_LB_LP_prep_for_ML.filter(regex='|'.join(featureSelection))
freeze5_LB_LP_prep_for_ML = freeze5_LB_LP_prep_for_ML.drop(featureRemoval, axis=1)

# Convert LB/LP labels into booleans
freeze5_LB_LP_prep_for_ML.loc[freeze5_LB_LP_prep_for_ML.ann_classificationVKGL == 'LB', 'ann_classificationVKGL'] = False
freeze5_LB_LP_prep_for_ML.loc[freeze5_LB_LP_prep_for_ML.ann_classificationVKGL == 'LP', 'ann_classificationVKGL'] = True
freeze5_LB_LP_prep_for_ML['ann_classificationVKGL'] = freeze5_LB_LP_prep_for_ML['ann_classificationVKGL'].astype('bool')
# Split into dependent (X) and independent (y) variables as well as train, validation and test sets
X = freeze5_LB_LP_prep_for_ML.drop('ann_classificationVKGL', axis=1)
y = freeze5_LB_LP_prep_for_ML['ann_classificationVKGL']
X_train, X_val_test, y_train, y_val_test = train_test_split(X, y, test_size=0.2, random_state=42)
X_val, X_test, y_val, y_test = train_test_split(X_val_test, y_val_test, test_size=0.5, random_state=42)


# Machine learning procedure
############################
# TODO
# implement 5 fold cross validation, perhaps stratification on y across folds
# Try different models (any or all of RF, GLMNET, Naive Bayes, Neaural Net, SVN, XGBoost, ElasticNet, SVM)
# Metrics: ROC AUC Mean, ROC AUC St.Dev, on tes/train, on 10-40-full features
# Find effect of regularization within models (e.g. in XGBoost L1, L2)
#
# Define XGBoost model, hyperparameters have been tuned on the validation test
xgb_reg = xgb.XGBRegressor(objective='binary:logistic', eval_metric = 'logloss', eta = 0.01, subsample = 0.6, learning_rate = 0.09, n_estimators = 150, random_state = 42)
# Using the train/validation data, find out optimal number of estimators
eval_set = [(X_train, y_train), (X_val, y_val)]
xgb_reg.fit(X_train, y_train, eval_set=eval_set)
# After optimizing, check the AUC of the validation set
print("XGBoost model validation AUC: ", roc_auc_score(y_val, xgb_reg.predict(X_val)))
# Incorporate the validation set into the fit to create the final model
xgb_reg.fit(pd.concat([X_train, X_val]), pd.concat([y_train, y_val]), eval_set=eval_set)
# After optimizing, check the AUC of the test set that is now seen for the first time
print("XGBoost model final test AUC: ", roc_auc_score(y_test, xgb_reg.predict(X_test)))


# Explain the predictions
#########################
# Set up SHAP explainer
explainer = shap.TreeExplainer(xgb_reg)
# sigmoid(SHAP) -> probability
# inverseSigmoid(probability) -> SHAP
# We set these anchor points:
# 1. base SHAP value (same for each observation) -> calculate probability
# 2. total SHAP value (different for each observation) -> calculate probability
# using these, we can use the per-feature SHAP value to scale probabilities
# NOTE:
# SHAP values are reported in the log-odds domain so that they preserve the additive property. Once you apply a sigmoid/logistic transform, you are in probability space—which is non-linear.
# While scaling each log-odds SHAP by a fraction of Δ𝑝 is a straightforward trick to “get some notion” of probability contribution, it does not correspond to the exact SHAP decomposition in probability space.
# It’s not wrong to say “I want to distribute the probability gap proportionally to each feature’s share of total log-odds SHAP.”
# But it is a heuristic (a “quick-and-dirty” approach) that often appears in practice, but it is not the strict, theory-backed approach that underpins SHAP.
# Functions:
def sigmoid(x):
    return 1 / (1 + np.exp(-x))
def inverseSigmoid(x: float):
    return np.log(x / (1 - x))

# Function to add scaled probabilities to input data
def predict_as_prob_scaled_SHAP(df: pd.DataFrame, explainer: shap.TreeExplainer, mlsplit: str) -> pd.DataFrame:
    expl = explainer(df)

    # Convert base SHAP value to probabilities
    df_SHAP_basevalue =  expl.base_values.reshape(-1, 1)
    df_SHAP_basevalue_P = sigmoid(df_SHAP_basevalue)

    # Sum per-feature SHAP values for each observations and convert to probabilities
    df_SHAP_values = expl.values
    df_SHAP_values_sum = df_SHAP_values.sum(axis=1, keepdims=True)
    df_SHAP_values_sum_P = sigmoid(df_SHAP_basevalue + df_SHAP_values_sum)

    # Sanity check: this should be equal to prediction probability output (except for rounding..)
    xgb_reg.predict(df).reshape(-1, 1)

    # Get scaling factors using SHAP values and apply to the probability range (= difference prob. sum and base value)
    df_SHAP_values_scaled = df_SHAP_values / df_SHAP_values_sum
    df_SHAP_values_scaled_P = df_SHAP_values_scaled * (df_SHAP_values_sum_P - df_SHAP_basevalue_P)

    # Rename and apply columns, add 'sph' for 'SHAP Probability Heuristic'
    renamed_columns = [f"{col}.sph" for col in df.columns]
    result = pd.DataFrame(df_SHAP_values_scaled_P, columns=renamed_columns, index=df.index)
    result["BaseProbability.sph"] = pd.DataFrame(df_SHAP_basevalue_P, index=df.index)
    result["FinalProbability.sph"] = pd.DataFrame(df_SHAP_values_sum_P, index=df.index)
    result["MLSplit"] = mlsplit
    return result

# Get predictions for X train, val, test as SHAP-value scales probabilities (rowsum = final probability)
X_train_SHAP_prob_values = predict_as_prob_scaled_SHAP(X_train, explainer, "Train")
X_val_SHAP_prob_values = predict_as_prob_scaled_SHAP(X_val, explainer, "Train")
X_test_SHAP_prob_values = predict_as_prob_scaled_SHAP(X_test, explainer, "Test")

# Prepare VUS (and CF) set for prediction and add probability scaled SHAP values
fr5_VUS_CF_forPred = freeze5_VUS_CF.filter(regex='|'.join(featureSelection))
fr5_VUS_CF_forPred = fr5_VUS_CF_forPred.drop(featureRemoval, axis=1)
fr5_VUS_CF_forPred = fr5_VUS_CF_forPred.drop('ann_classificationVKGL', axis=1)
fr5_VUS_CF_SHAP_prob_values = predict_as_prob_scaled_SHAP(fr5_VUS_CF_forPred, explainer, "Unknown")

# Combine all with original data and write to file
allPredsInOne = pd.concat([X_train_SHAP_prob_values, X_val_SHAP_prob_values, X_test_SHAP_prob_values, fr5_VUS_CF_SHAP_prob_values])
freeze5_plus_pred = pd.concat([freeze5, allPredsInOne], axis=1)
freeze5_plus_pred.to_csv("/Users/joeri/git/dave/data/freeze5_predictions.csv", index=False)

