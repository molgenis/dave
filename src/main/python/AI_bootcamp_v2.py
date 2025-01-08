# Use Python 3.12.4
import pandas as pd
import numpy as np
import xgboost as xgb
import matplotlib.pyplot as plt
import graphviz
import shap
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score, confusion_matrix, precision_score, recall_score, ConfusionMatrixDisplay, precision_recall_curve
from sklearn.model_selection import RandomizedSearchCV, train_test_split, GridSearchCV
from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor, plot_tree
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from scipy.stats import randint
from sklearn.tree import export_graphviz
from sklearn import tree
from IPython.display import Image
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC
from IPython.core.display import display, HTML

# Solving this error: The XGBRegressor or classes from which it inherits use `_get_tags` and `_more_tags`
# https://stackoverflow.com/questions/79290968/super-object-has-no-attribute-sklearn-tags
# --> downgrade scikit-learn
# pip uninstall -y scikit-learn
# pip install scikit-learn==1.5.2

# Load source data and select only columns that are good for ML
# add 'ann_am_pathogenicity' to train upon AlphaMissense
# convert protein localization into booleans
freeze5 = pd.read_csv("/Users/joeri/git/vkgl-secretome-protein-stability/data/freeze5-provisional.csv.gz", low_memory=False)
selectColumns = ["ann_classificationVKGL", "ann_proteinIschaperoned", "ann_proteinLocalization", "delta_"]
freeze5 = freeze5.filter(regex='|'.join(selectColumns))
freeze5 = freeze5.drop(['delta_aaSeq'], axis=1)
freeze5 = pd.get_dummies(freeze5, columns=['ann_proteinLocalization'], drop_first=True)

# For modeling purposes, create 'fr5_ML' with only binarized LP/LB mutations, and dummies for protein localization
fr5_ML = freeze5[(freeze5['ann_classificationVKGL'] == 'LP') | (freeze5['ann_classificationVKGL'] == 'LB')]
fr5_ML.loc[fr5_ML.ann_classificationVKGL == 'LB', 'ann_classificationVKGL'] = False
fr5_ML.loc[fr5_ML.ann_classificationVKGL == 'LP', 'ann_classificationVKGL'] = True
fr5_ML['ann_classificationVKGL'] = fr5_ML['ann_classificationVKGL'].astype('bool')
#fr5_ML.info(verbose=True)

######
# ML #
######

# Split in dep/indep vars and test/train sets, train RF and determine AUC
X = fr5_ML.drop('ann_classificationVKGL', axis=1)
y = fr5_ML['ann_classificationVKGL']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# XGBoost
eval_set = [(X_train, y_train), (X_test, y_test)]
xgb_reg = xgb.XGBRegressor(objective='binary:logistic', eval_metric = 'logloss', eta = 0.01, subsample = 0.6, learning_rate = 0.09, n_estimators = 150, random_state = 42)
xgb_reg.fit(X_train, y_train, eval_set=eval_set)
print("XGBoost model AUC: ", roc_auc_score(y_test, xgb_reg.predict(X_test))) # 0.9004719623873337

# Plot fit per boosting round
results = xgb_reg.evals_result()
epochs = len(results['validation_0']['logloss'])
x_axis = range(0, epochs)
plt.figure(figsize=(10, 6))
plt.plot(x_axis, results['validation_0']['logloss'], label='Train')
plt.plot(x_axis, results['validation_1']['logloss'], label='Test')
plt.xlabel('Boosting Rounds')
plt.ylabel('Log Loss')
plt.title('XGBoost Log Loss')
plt.legend()
plt.grid(True)
plt.show()

# Plot AUC
fpr, tpr, thresholds = roc_curve(y_test, xgb_reg.predict(X_test))
plt.figure()
plt.plot(fpr, tpr, color='blue', lw=2, label=f'ROC curve (AUC = {roc_auc_score(y_test, xgb_reg.predict(X_test)):.2f})')
plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate (FPR)')
plt.ylabel('True Positive Rate (TPR)')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc="lower right")
plt.show()

# Feature importance
featureImpResult = permutation_importance(xgb_reg, X_test, y_test, n_repeats=10, random_state=42, n_jobs=4)
forest_importances = pd.Series(featureImpResult.importances_mean, index=X_train.columns)
forest_importances.sort_values(ascending=False)

# SHAP beeswarm of feature impact
explainer = shap.TreeExplainer(xgb_reg)
shap_values = explainer.shap_values(X)
shap.summary_plot(shap_values, X)

############################
# Apply to VUS and explain #
############################

# Add probabilities to all original rows? niet handig ivm bewerkte feature names?
#probAll = xgb_reg.predict(freeze5.drop("ann_classificationVKGL", axis=1))
#freeze5["Prediction_probability"] = probAll

# Function to calculate PPV/NPV for a given threshold
def calculate_metric(y_true, y_pred, threshold, metric):
    preds = y_pred >= threshold
    tp = np.sum((preds == 1) & (y_true == 1))
    fp = np.sum((preds == 1) & (y_true == 0))
    tn = np.sum((preds == 0) & (y_true == 0))
    fn = np.sum((preds == 0) & (y_true == 1))
    if metric == "ppv":
        return tp / (tp + fp) if (tp + fp) > 0 else 0
    elif metric == "npv":
        return tn / (tn + fn) if (tn + fn) > 0 else 0
    else:
        raise ValueError("Unsupported metric. Use 'ppv' or 'npv'.")

# Add probabilities to test set and find PPV/NPV thresholds
probs = xgb_reg.predict(X_test)
thresholds = np.linspace(0, 1, 100)
ppv_values = [calculate_metric(y_test, probs, t, metric="ppv") for t in thresholds]
npv_values = [calculate_metric(y_test, probs, t, metric="npv") for t in thresholds]

# Find threshold with >= 95% PPV
threshold =  0.90
threshold_ppv = thresholds[np.where(np.array(ppv_values) >= threshold)[0][0]]
threshold_npv = thresholds[np.where(np.array(npv_values) >= threshold)[0][0]]
print("Threshold for",threshold*100,"% PPV:", threshold_ppv)
print("Threshold for",threshold*100,"% NPV:", threshold_npv)

# Predict probabilities for VUS observations
fr5_VUS = freeze5[(freeze5['ann_classificationVKGL'] == 'VUS')]
fr5_VUS = fr5_VUS[X_train.columns]  # quick way to align columns (e.g. remove "ann_classificationVKGL")
vus_probs = xgb_reg.predict(fr5_VUS)

# Filter VUS with predictions above the PPV and below NPV threshold
confident_LP = fr5_VUS[vus_probs >= threshold_ppv]
confident_LP['Predicted_Label'] = "LP"
confident_LP['Prediction_Probability'] = vus_probs[vus_probs >= threshold_ppv]
confident_LP.shape

confident_LB = fr5_VUS[vus_probs <= threshold_npv]
confident_LB['Predicted_Label'] = "LB"
confident_LB['Prediction_Probability'] = vus_probs[vus_probs <= threshold_npv]
confident_LB.shape

confidentPred = pd.concat([confident_LB, confident_LP], axis=0)
confidentPred.shape

# SHAP explanation for each confident VUS observation
# see: https://shap.readthedocs.io/en/latest/tabular_examples.html
explainer = shap.TreeExplainer(xgb_reg)
shapX = confidentPred.drop(['Predicted_Label', 'Prediction_Probability'], axis=1)
shap_values = explainer.shap_values(shapX)

# Display SHAP explanations
i = 1
one_force_plot = shap.force_plot(explainer.expected_value, shap_values[i], confidentPred.drop(['Predicted_Label', 'Prediction_Probability'], axis=1).iloc[i], matplotlib=False)
one_shap_html = f"{shap.getjs()}{one_force_plot.html()}"
display(HTML(one_shap_html))
file = open("one_shapley.html", "w")
file.write(one_shap_html)

all_force_plot = shap.force_plot(explainer.expected_value, shap_values[:shap_values.shape[0], :], confidentPred.drop(['Predicted_Label', 'Prediction_Probability'], axis=1).iloc[:shap_values.shape[0], :], matplotlib=False)
all_shap_html = f"{shap.getjs()}{all_force_plot.html()}"
display(HTML(all_shap_html))
file = open("shapley.html", "w")
file.write(all_shap_html)

# Save results
#confidentPred.to_csv("confident_vus_predictions.csv", index=False)
