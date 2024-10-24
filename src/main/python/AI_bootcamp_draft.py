# Data Processing
import pandas as pd
import numpy as np

# Modelling
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score, ConfusionMatrixDisplay
from sklearn.model_selection import RandomizedSearchCV, train_test_split
from scipy.stats import randint

# Tree Visualisation
from sklearn.tree import export_graphviz
from IPython.display import Image
import graphviz
from sklearn.metrics import roc_auc_score

dat = pd.read_csv("data/freeze4.csv.gz", low_memory=False)

# "ann_proteinLocalization",
# small model
#search_strings = ["ann_classificationVKGL", "ann_proteinIschaperoned", "delta_aliphaticIndex", "delta_backbone.clash", "delta_Backbone.Hbond", "delta_bomanIndex", "delta_charge", "delta_cis_bond", "delta_disulfide","delta_Electrostatics", "delta_energy.Ionisation", "delta_entropy.mainchain", "delta_entropy.sidechain", "delta_helix.dipole", "delta_hydrophobicity", "delta_hydrophobicMoment", "delta_instabilityIndex", "delta_isoElecPoint", "delta_massOverCharge", "delta_molWeight", "delta_Sidechain.Hbond", "delta_Solvation.Hydrophobic", "delta_Solvation.Polar", "delta_torsional.clash", "delta_total.energy", "delta_Van.der.Waals", "delta_Van.der.Waals.clashes"]
# big model
search_strings = ["ann_classificationVKGL", "ann_proteinIschaperoned", "ann_proteinLocalization", "delta_"]
dat = dat.filter(regex='|'.join(search_strings))

# Drop columns
dat.drop(columns=['delta_aaSeq'], inplace=True)

# Filter on LB/LP and binarize target variable
dat = dat[(dat['ann_classificationVKGL'] == 'LP') | (dat['ann_classificationVKGL'] == 'LB')]
dat['ann_classificationVKGL'] = dat['ann_classificationVKGL'].replace({'LB': 0, 'LP': 1})

# Add dummy variabless
dat = pd.get_dummies(dat, columns=['ann_proteinLocalization'], drop_first=True)

# data is now ready for analysis!!
print("columns = ", dat.columns)
print("shape = ", dat.shape)

# Split the data into features (X) and target (y)
X = dat.drop('ann_classificationVKGL', axis=1)
y = dat['ann_classificationVKGL']
print("X = ", X.shape, " y = ", y.shape)

# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

# fit RF
rf = RandomForestClassifier(
    n_estimators = 100,     # 100 - 500 --> MAKES NO DIFFERENCE
    min_samples_split = 10,  # 2 - 10 --> MAKES NO DIFFERENCE
    min_samples_leaf = 6,   # 1 - 6 --> MAKES NO DIFFERENCE
    max_features = 'sqrt',  # sqrt/log2 --> sqrt slightly better
    max_depth = 25,       # None or 10 - 25 --> higher slightly better
    random_state = 42       # Random seed for reproducibility
)
rf.fit(X_train, y_train)

# probabilities of unknowns (Probability for class 1)
y_probs = rf.predict_proba(X_test)[:, 1]

# predict classes of the unknowns
#y_pred = rf.predict(X_test)

# AUC
from sklearn.metrics import roc_auc_score
auc_score = roc_auc_score(y_test, y_probs)
print(f"AUC: {auc_score:.3f}")


# plot AUC
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.datasets import make_classification
# 5. Compute and plot the ROC curve
fpr, tpr, thresholds = roc_curve(y_test, y_probs)
plt.figure()
plt.plot(fpr, tpr, color='blue', lw=2, label=f'ROC curve (AUC = {auc_score:.2f})')
plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')  # Dashed diagonal
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate (FPR)')
plt.ylabel('True Positive Rate (TPR)')
plt.title('Receiver Operating Characteristic (ROC) Curve')
plt.legend(loc="lower right")
plt.show()

# permuted feat imp --> use these to priorize features in next iteration
from sklearn.inspection import permutation_importance
result = permutation_importance(rf, X_test, y_test, n_repeats=10, random_state=42, n_jobs=2)
forest_importances = pd.Series(result.importances_mean, index=X_train.columns)
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots()
# forest_importances.plot.bar(yerr=result.importances_std, ax=ax)
# ax.set_title("Feature importances using permutation on full model")
# ax.set_ylabel("Mean accuracy decrease")
# fig.tight_layout()
# plt.show()

# print top X features
topX = 25
top_features = forest_importances.sort_values(ascending=False).head(topX)
print(top_features)

# prepare top X feature list and filter

# filter on top features, train next model
topFeatDat = dat[['ann_classificationVKGL'] + top_features.index.tolist()]
X = topFeatDat.drop('ann_classificationVKGL', axis=1)
y = topFeatDat['ann_classificationVKGL']
print("X = ", X.shape, " y = ", y.shape)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

rf = RandomForestClassifier(
    n_estimators = 200,      # 100 - 500
    min_samples_split = 6,   # 2 - 10 --> MAKES NO DIFFERENCE
    min_samples_leaf = 2,    # 1 - 6  --> MAKES NO DIFFERENCE
    max_features = 'log2',   # sqrt/log2 --> log2 slighly better
    max_depth = 25,          # None or 10 - 25 --> higher is slighly better
    random_state = 42        # Random seed for reproducibility
)
rf.fit(X_train, y_train)
y_probs = rf.predict_proba(X_test)[:, 1]
from sklearn.metrics import roc_auc_score
auc_score = roc_auc_score(y_test, y_probs)
print(f"AUC: {auc_score:.3f}")


##### XGBoost!!

import xgboost as xgb

topFeatDat = dat[['ann_classificationVKGL'] + top_features.index.tolist()]
X = topFeatDat.drop('ann_classificationVKGL', axis=1)
y = topFeatDat['ann_classificationVKGL']

print("X = ", X.shape, " y = ", y.shape)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
print(y_train)


# data_dmatrix = xgb.DMatrix(data=X_train,label=y_train)
# xgb_cv = xgb.cv(dtrain=data_dmatrix, params=params, nfold=5, metrics = 'logloss',seed=42) 

xgb_reg = xgb.XGBRegressor(objective='binary:logistic',
                           eval_metric = 'logloss',
                           eta = 0.01, # no diff
                           subsample = 0.1, # no diff
                           learning_rate = 0.1, # no diff
                           n_estimators = 10, # no diff
                           use_label_encoder = True) # no diff
xgb_reg.fit(X_train, y_train) # nfold ?


# Make predictions on the test set
y_pred_proba = xgb_reg.predict(X_test)

# Calculate the ROC AUC score
roc_auc = roc_auc_score(y_test, y_pred_proba)
print(f"ROC AUC Score: {roc_auc:.2f}")

import shap
explainer = shap.TreeExplainer(xgb_reg)
explanation = explainer(X_test)
shap_values = explanation.values
shap.plots.beeswarm(explanation)
shap.initjs()
import matplotlib
shap.force_plot(explainer.expected_value, shap_values[0,:], X_test.iloc[0,:], matplotlib=matplotlib)


### DECISION TREE

from sklearn.tree import DecisionTreeClassifier # Import Decision Tree Classifier
topFeatDat = dat[['ann_classificationVKGL'] + top_features.index.tolist()]
X = topFeatDat.drop('ann_classificationVKGL', axis=1)
y = topFeatDat['ann_classificationVKGL']
print("X = ", X.shape, " y = ", y.shape)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Create Decision Tree classifer object
clf = DecisionTreeClassifier(criterion="gini", max_depth=6)

# Train Decision Tree Classifer
clf = clf.fit(X_train,y_train)

#Predict the response for test dataset
y_pred_proba = clf.predict_proba(X_test)[:, 1]

roc_auc = roc_auc_score(y_test, y_pred_proba)
print(f"ROC AUC Score: {roc_auc:.2f}")



from sklearn.tree import export_graphviz
from sklearn.externals.six import StringIO  
from IPython.display import Image  
import pydotplus

dot_data = StringIO()
export_graphviz(clf, out_file=dot_data,  
                filled=True, rounded=True,
                special_characters=True,feature_names = X_test.columns,class_names=['0','1'])
graph = pydotplus.graph_from_dot_data(dot_data.getvalue())  
graph.write_png('diabetes.png')
Image(graph.create_png())