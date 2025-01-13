# Use Python 3.12.4
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score


# Import data
freeze5 = pd.read_csv("/Users/joeri/git/vkgl-secretome-protein-stability/data/freeze5-provisional-unbiased-impute.csv.gz", low_memory=False)
# Split once into LB/LP and VUS/CF and do not touch
freeze5_LB_LP = freeze5[(freeze5['ann_classificationVKGL'] == 'LB') | (freeze5['ann_classificationVKGL'] == 'LP')]
freeze5_VUS_CF = freeze5[(freeze5['ann_classificationVKGL'] == 'VUS') | (freeze5['ann_classificationVKGL'] == 'CF')]
# i.e. backTogether = pd.concat([freeze5_LB_LP, freeze5_VUS_CF], axis=0)

# Data preparation: keep indices of LB/LP set the same so we can merge results later
# Change protein localization from 3 categoricals into 2 booleans: "xx_membrane" and
# "xx_secreted" with the 3rd one ("xx_intracellular") being implicit if the other ones are false
freeze5_LB_LP_prep_for_ML = pd.get_dummies(freeze5_LB_LP, columns=['ann_proteinLocalization'], drop_first=True)
# Select columns that are relevant features for training
featureSelection = ["ann_classificationVKGL", "ann_proteinIschaperoned", "ann_proteinLocalization", "delta_"]
freeze5_LB_LP_prep_for_ML = freeze5_LB_LP_prep_for_ML.filter(regex='|'.join(featureSelection))
freeze5_LB_LP_prep_for_ML = freeze5_LB_LP_prep_for_ML.drop(['delta_aaSeq'], axis=1)
# For modeling purposes, make booleans for of LB andLP labels
freeze5_LB_LP_prep_for_ML.loc[freeze5_LB_LP_prep_for_ML.ann_classificationVKGL == 'LB', 'ann_classificationVKGL'] = False
freeze5_LB_LP_prep_for_ML.loc[freeze5_LB_LP_prep_for_ML.ann_classificationVKGL == 'LP', 'ann_classificationVKGL'] = True
freeze5_LB_LP_prep_for_ML['ann_classificationVKGL'] = freeze5_LB_LP_prep_for_ML['ann_classificationVKGL'].astype('bool')

# Split in dep/indep vars and test/train sets, train RF and determine AUC
X = freeze5_LB_LP_prep_for_ML.drop('ann_classificationVKGL', axis=1)
y = freeze5_LB_LP_prep_for_ML['ann_classificationVKGL']
X_train, X_val_test, y_train, y_val_test = train_test_split(X, y, test_size=0.2, random_state=42)
X_val, X_test, y_val, y_test = train_test_split(X_val_test, y_val_test, test_size=0.5, random_state=42)

### hyperparam tuning / Kfold CV ?
eval_set = [(X_train, y_train), (X_val, y_val)]
# Optimize params and boosting rounds using validation set
xgb_reg = xgb.XGBRegressor(objective='binary:logistic', eval_metric = 'logloss', eta = 0.01, subsample = 0.6, learning_rate = 0.09, n_estimators = 150, random_state = 42)
xgb_reg.fit(X_train, y_train, eval_set=eval_set)
print("XGBoost model validation AUC: ", roc_auc_score(y_val, xgb_reg.predict(X_val)))
# Incorporate validation set into the fit for final model
xgb_reg.fit(pd.concat([X_train, X_val]), pd.concat([y_train, y_val]), eval_set=eval_set)
print("XGBoost model final test AUC: ", roc_auc_score(y_test, xgb_reg.predict(X_test)))

# get probs back in original, how?
X_train["Probs"] = xgb_reg.predict(X_train)
X_train["MLSplit"] = "Train"
X_val["Probs"] = xgb_reg.predict(X_val)
X_val["MLSplit"] = "Validation"
X_test["Probs"] = xgb_reg.predict(X_test)
X_test["MLSplit"] = "Test"
X_all = pd.concat([X_train, X_val, X_test])
#X_all.sort_index(axis=0, inplace=True) #original X back, sort of

freeze5_VUS_CF_forPred = pd.get_dummies(freeze5_VUS_CF, columns=['ann_proteinLocalization'], drop_first=True)
freeze5_VUS_CF_forPred = freeze5_VUS_CF_forPred.filter(regex='|'.join(featureSelection))
freeze5_VUS_CF_forPred = freeze5_VUS_CF_forPred.drop(['delta_aaSeq'], axis=1)
freeze5_VUS_CF_forPred = freeze5_VUS_CF_forPred.drop('ann_classificationVKGL', axis=1)
freeze5_VUS_CF_forPred["Probs"] = xgb_reg.predict(freeze5_VUS_CF_forPred)
freeze5_VUS_CF_forPred["MLSplit"] = "Target"

allPredsInOne = pd.concat([X_all, freeze5_VUS_CF_forPred])
freeze5["Probs"] = allPredsInOne["Probs"]
freeze5["MLSplit"] = allPredsInOne["MLSplit"]
freeze5.to_csv("/Users/joeri/git/vkgl-secretome-protein-stability/data/freeze5_with_predictions.csv", index=False)
