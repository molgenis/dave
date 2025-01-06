import pandas as pd
import numpy as np
import xgboost as xgb
import matplotlib.pyplot as plt
import graphviz
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score, confusion_matrix, precision_score, recall_score, ConfusionMatrixDisplay
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
#pd.set_option('future.no_silent_downcasting', True)

# Solving this error: The XGBRegressor or classes from which it inherits use `_get_tags` and `_more_tags`
# https://stackoverflow.com/questions/79290968/super-object-has-no-attribute-sklearn-tags
# --> downgrade scikit-learn
# pip uninstall -y scikit-learn
# pip install scikit-learn==1.5.2

# Load source data and select only columns that are good for ML
# add 'ann_am_pathogenicity' to train upon AlphaMissense
freeze5 = pd.read_csv("/Users/joeri/git/vkgl-secretome-protein-stability/data/freeze5-provisional.csv.gz", low_memory=False)
anyInterestingColumns = ["ann_classificationVKGL", "ann_proteinIschaperoned", "ann_proteinLocalization", "delta_"]
freeze5 = freeze5.filter(regex='|'.join(anyInterestingColumns))
freeze5 = freeze5.drop(['delta_aaSeq'], axis=1)

# For modeling purposes, create 'fr5_ML' with only binarized LP/LB mutations, and dummies for protein localization
fr5_ML = freeze5[(freeze5['ann_classificationVKGL'] == 'LP') | (freeze5['ann_classificationVKGL'] == 'LB')]
fr5_ML.loc[fr5_ML.ann_classificationVKGL == 'LB', 'ann_classificationVKGL'] = False
fr5_ML.loc[fr5_ML.ann_classificationVKGL == 'LP', 'ann_classificationVKGL'] = True
fr5_ML['ann_classificationVKGL'] = fr5_ML['ann_classificationVKGL'].astype('bool')
fr5_ML = pd.get_dummies(fr5_ML, columns=['ann_proteinLocalization'], drop_first=True)
#fr5_ML.info(verbose=True)

######
# ML #
######

# Split in dep/indep vars and test/train sets, train RF and determine AUC
X = fr5_ML.drop('ann_classificationVKGL', axis=1)
y = fr5_ML['ann_classificationVKGL']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Random Forest, always a good baseline
#rf = RandomForestClassifier(n_estimators = 150, min_samples_split = 6, min_samples_leaf = 2, max_features = 'sqrt', max_depth = 20, random_state = 42)
#rf.fit(X_train, y_train)
#print("Random forest full feature model AUC:", roc_auc_score(y_test, rf.predict_proba(X_test)[:, 1])) # 0.87656

# XGBoost
eval_set = [(X_train, y_train), (X_test, y_test)]
xgb_reg = xgb.XGBRegressor(objective='binary:logistic', eval_metric = 'logloss', eta = 0.01, subsample = 0.6, learning_rate = 0.09, n_estimators = 150, random_state = 42)
xgb_reg.fit(X_train, y_train, eval_set=eval_set)
print("XGBoost model AUC: ", roc_auc_score(y_test, xgb_reg.predict(X_test))) # 0.9004719623873337

# Plot fit per boosting round #
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

# Plot AUC #
fpr, tpr, thresholds = roc_curve(y_test, xgb_reg.predict(X_test))
plt.figure()
plt.plot(fpr, tpr, color='blue', lw=2, label=f'ROC curve (AUC = {roc_auc_score(y_test, xgb_reg.predict(X_test)):.2f})')
plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')  # Dashed diagonal
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate (FPR)')
plt.ylabel('True Positive Rate (TPR)')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc="lower right")
plt.show()