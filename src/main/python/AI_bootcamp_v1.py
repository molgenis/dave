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


# Load source data
df_source = pd.read_csv("data/freeze4.csv.gz", low_memory=False)

# Select only columns that are potentially interesting and replace df_source to save memory
anyInterestingColumns = ["ann_classificationVKGL", "ann_proteinIschaperoned", "ann_am_pathogenicity", "ann_proteinLocalization", "delta_"]
df_source = df_source.filter(regex='|'.join(anyInterestingColumns))
df_source.drop(columns=['delta_aaSeq'], inplace=True)

# For modeling purposes, create 'df_forML' with only binarized LP/LB mutations, and dummies for localization
df_forML = df_source[(df_source['ann_classificationVKGL'] == 'LP') | (dat['ann_classificationVKGL'] == 'LB')]
df_forML['ann_classificationVKGL'] = df_forML['ann_classificationVKGL'].replace({'LB': 0, 'LP': 1})
df_forML = pd.get_dummies(df_forML, columns=['ann_proteinLocalization'], drop_first=True)
print("columns = ", df_forML.columns)
print("shape = ", df_forML.shape)
print(df_forML.columns)

# Optional data selection: e.g. only secreted proteins
#df_forML = df_forML[(df_forML['ann_proteinLocalization_secreted'] == 1)]

###################################
# Question 1: Can we outperform AM?
###################################

# Drop AlphaMissense for this question
df_noAM = df_forML.drop(columns=['ann_am_pathogenicity'], inplace=False)

# Split in dep/indep vars and test/train sets, train RF and determine AUC
X = df_noAM.drop('ann_classificationVKGL', axis=1)
y = df_noAM['ann_classificationVKGL']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Random Forest, always a good baseline
rf = RandomForestClassifier(n_estimators = 150, min_samples_split = 6, min_samples_leaf = 2, max_features = 'sqrt', max_depth = 20, random_state = 42)
rf.fit(X_train, y_train)
print("Random forest full feature model AUC:", roc_auc_score(y_test, rf.predict_proba(X_test)[:, 1])) # ~ .86

# XGBoost
xgb_reg = xgb.XGBRegressor(objective='binary:logistic', eval_metric = 'logloss', eta = 0.01, subsample = 0.6, learning_rate = 0.1, n_estimators = 75, random_state = 42)
xgb_reg.fit(X_train, y_train)
print("XGBoost full feature model AUC: ", roc_auc_score(y_test, xgb_reg.predict(X_test))) # ~ .88 (87.5392)

# Logistic regression
logreg = LogisticRegression()
logreg.fit(X_train, y_train)
print("Logistic regression full feature model AUC", roc_auc_score(y_test, logreg.predict_proba(X_test)[:, 1])) # ~ .78

# Multi-layer Perceptron
mlp_cl = MLPClassifier()
mlp_cl.fit(X_train, y_train)
print("Multi-layer Perceptron full feature model AUC", roc_auc_score(y_test, mlp_cl.predict_proba(X_test)[:, 1])) # ~ .85

# Support Vector Machine
svm_cl = SVC(probability=True)
svm_cl.fit(X_train, y_train)
print("Support Vector Machine full feature model AUC", roc_auc_score(y_test, svm_cl.predict_proba(X_test)[:, 1])) # ~ .79

#######

# From the RF model, determine feature importance by permutation and select top features
featureImpResult = permutation_importance(rf, X_test, y_test, n_repeats=10, random_state=42, n_jobs=2)
forest_importances = pd.Series(featureImpResult.importances_mean, index=X_train.columns)
top_features = forest_importances.sort_values(ascending=False).head(25)

# Re-train RF on top features only
df_noAM_topF = df_noAM[['ann_classificationVKGL'] + top_features.index.tolist()]
X = df_noAM_topF.drop('ann_classificationVKGL', axis=1)
y = df_noAM_topF['ann_classificationVKGL']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
rf = RandomForestClassifier(n_estimators = 150, min_samples_split = 2, min_samples_leaf = 2, max_features = 'sqrt', max_depth = 15, random_state = 42)
rf.fit(X_train, y_train)
print("Random forest top feature model AUC:", roc_auc_score(y_test, rf.predict_proba(X_test)[:, 1])) # ~ .84

# Try XGBoost on the same set
xgb_reg = xgb.XGBRegressor(objective='binary:logistic', eval_metric = 'logloss', eta = 0.01, subsample = 0.5, learning_rate = 0.1, n_estimators = 75, random_state = 42)
xgb_reg.fit(X_train, y_train)
print("XGBoost top feature model AUC: ", roc_auc_score(y_test, xgb_reg.predict(X_test))) # ~ .84

# Try logistic regression
logreg = LogisticRegression()
logreg.fit(X_train, y_train)
print("Logistic regression top feature model AUC:", roc_auc_score(y_test, logreg.predict_proba(X_test)[:, 1])) # ~ .72

# Multi-layer Perceptron
mlp_cl = MLPClassifier()
mlp_cl.fit(X_train, y_train)
print("Multi-layer Perceptron top feature model AUC", roc_auc_score(y_test, mlp_cl.predict_proba(X_test)[:, 1])) # ~ .80

# Support Vector Machine
svm_cl = SVC(probability=True)
svm_cl.fit(X_train, y_train)
print("Support Vector Machine top feature model AUC", roc_auc_score(y_test, svm_cl.predict_proba(X_test)[:, 1])) # ~ .72

####################################################
# Question 2: Can we combine with AM and do better?
####################################################

# drop rows without ann_am_pathogenicity
df_rowsWithAM = df_forML.dropna(subset=['ann_am_pathogenicity'], inplace=False)

# Split in dep/indep vars and test/train sets, train RF and determine AUC
X = df_rowsWithAM.drop('ann_classificationVKGL', axis=1)
y = df_rowsWithAM['ann_classificationVKGL']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
rf = RandomForestClassifier(n_estimators = 150, min_samples_split = 6, min_samples_leaf = 2, max_features = 'sqrt', max_depth = 20, random_state = 42)
rf.fit(X_train, y_train)
print("Random forest full feature model AUC:", roc_auc_score(y_test, rf.predict_proba(X_test)[:, 1])) # 0.93637

# XGBoost
xgb_reg = xgb.XGBRegressor(objective='binary:logistic', eval_metric = 'logloss', eta = 0.01, subsample = 0.6, learning_rate = 0.1, n_estimators = 75, random_state = 42)
xgb_reg.fit(X_train, y_train)
print("XGBoost full feature model AUC: ", roc_auc_score(y_test, xgb_reg.predict(X_test))) #0.955505

# Logistic regression
logreg = LogisticRegression()
logreg.fit(X_train, y_train)
print("Logistic regression full feature model AUC", roc_auc_score(y_test, logreg.predict_proba(X_test)[:, 1])) # 0.94225

# Multi-layer Perceptron
mlp_cl = MLPClassifier()
mlp_cl.fit(X_train, y_train)
print("Multi-layer Perceptron top feature model AUC", roc_auc_score(y_test, mlp_cl.predict_proba(X_test)[:, 1])) # 0.94782

# Support Vector Machine
svm_cl = SVC(probability=True)
svm_cl.fit(X_train, y_train)
print("Support Vector Machine top feature model AUC", roc_auc_score(y_test, svm_cl.predict_proba(X_test)[:, 1])) # 0.88166

# From the RF model, determine feature importance by permutation and select top features
featureImpResult = permutation_importance(rf, X_test, y_test, n_repeats=10, random_state=42, n_jobs=2)
forest_importances = pd.Series(featureImpResult.importances_mean, index=X_train.columns)
top_features = forest_importances.sort_values(ascending=False).head(25)

# Re-train RF on top features only
df_rowsWithAM_topF = df_rowsWithAM[['ann_classificationVKGL'] + top_features.index.tolist()]
X = df_rowsWithAM_topF.drop('ann_classificationVKGL', axis=1)
y = df_rowsWithAM_topF['ann_classificationVKGL']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
rf = RandomForestClassifier(n_estimators = 150, min_samples_split = 2, min_samples_leaf = 2, max_features = 'sqrt', max_depth = 15, random_state = 42)
rf.fit(X_train, y_train)
print("Random forest top feature model AUC:", roc_auc_score(y_test, rf.predict_proba(X_test)[:, 1])) # 0.94578

# XGBoost on the same set
xgb_reg = xgb.XGBRegressor(objective='binary:logistic', eval_metric = 'logloss', eta = 0.01, subsample = 0.6, learning_rate = 0.1, n_estimators = 75, random_state = 42)
xgb_reg.fit(X_train, y_train)
print("XGBoost top feature model AUC: ", roc_auc_score(y_test, xgb_reg.predict(X_test))) # 0.952093

# Logistic regression
logreg = LogisticRegression()
logreg.fit(X_train, y_train)
print("Logistic regression top feature model AUC:", roc_auc_score(y_test, logreg.predict_proba(X_test)[:, 1])) # 0.94206

# Multi-layer Perceptron
mlp_cl = MLPClassifier()
mlp_cl.fit(X_train, y_train)
print("Multi-layer Perceptron top feature model AUC", roc_auc_score(y_test, mlp_cl.predict_proba(X_test)[:, 1])) # 0.94765

# Support Vector Machine
svm_cl = SVC(probability=True)
svm_cl.fit(X_train, y_train)
print("Support Vector Machine top feature model AUC", roc_auc_score(y_test, svm_cl.predict_proba(X_test)[:, 1])) # 0.93627

####################################################
# Question 3: Can we develop a viable predictor that is more insightful than AlphaMissense?
####################################################
df_explF = df_forML[['ann_classificationVKGL', 'ann_proteinIschaperoned', 'ann_proteinLocalization_membrane', 'ann_proteinLocalization_secreted', 'delta_aliphaticIndex', 'delta_backbone.clash', 'delta_Backbone.Hbond', 'delta_bomanIndex', 'delta_charge', 'delta_cis_bond', 'delta_disulfide', 'delta_electrostatic.kon', 'delta_Electrostatics', 'delta_energy.Ionisation', 'delta_Entropy.Complex', 'delta_entropy.mainchain', 'delta_entropy.sidechain', 'delta_helix.dipole', 'delta_hydrophobicity', 'delta_hydrophobicMoment', 'delta_instabilityIndex', 'delta_isoElecPoint', 'delta_massOverCharge', 'delta_mloop_entropy', 'delta_molWeight', 'delta_partial.covalent.bonds', 'delta_Sidechain.Hbond', 'delta_sloop_entropy', 'delta_Solvation.Hydrophobic', 'delta_Solvation.Polar', 'delta_torsional.clash', 'delta_total.energy', 'delta_Van.der.Waals', 'delta_Van.der.Waals.clashes', 'delta_water.bridge']]
X = df_explF.drop('ann_classificationVKGL', axis=1)
y = df_explF['ann_classificationVKGL']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
dt = DecisionTreeClassifier(random_state=42)
# Define the parameter grid for GridSearchCV
param_grid = {
    'max_depth': [3, 5, 10, None],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4],
    'max_features': ['sqrt', 'log2', None],
    'criterion': ['gini', 'entropy']
}
grid_search = GridSearchCV(estimator=dt, param_grid=param_grid, cv=5, verbose=1, n_jobs=-1)
grid_search.fit(X_train, y_train)
print("Best parameters found: ", grid_search.best_params_)
dt = DecisionTreeClassifier(random_state=42, criterion='gini', max_depth=5, max_features=None, min_samples_leaf=1, min_samples_split=2)
dt.fit(X_train, y_train)
print("Decision tree with explainable feature model AUC:", roc_auc_score(y_test, dt.predict_proba(X_test)[:, 1])) # 0.76270635

# Decision tree plot
dot_data = tree.export_graphviz(dt, out_file=None, feature_names=X.columns, class_names=['Pathogenic', 'Benign'], filled=True)
graph = graphviz.Source(dot_data, format="png")
graph.render("decision_tree_graph", format="png", cleanup=True)


######################################
# Check overfit on...
######################################

### either the 'full XGB model'
df_rowsWithAM = df_forML.dropna(subset=['ann_am_pathogenicity'], inplace=False)
X = df_rowsWithAM.drop('ann_classificationVKGL', axis=1)
y = df_rowsWithAM['ann_classificationVKGL']

### or top feature model
df_rowsWithAM_topF = df_rowsWithAM[['ann_classificationVKGL'] + top_features.index.tolist()]
X = df_rowsWithAM_topF.drop('ann_classificationVKGL', axis=1)
y = df_rowsWithAM_topF['ann_classificationVKGL']

### followed by
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
eval_set = [(X_train, y_train), (X_test, y_test)]
xgb_reg = xgb.XGBRegressor(objective='binary:logistic', eval_metric = 'logloss', eta = 0.01, subsample = 0.6, learning_rate = 0.1, n_estimators = 75, random_state = 42)
xgb_reg.fit(X_train, y_train, eval_set=eval_set, verbose=True)
print("XGBoost model AUC: ", roc_auc_score(y_test, xgb_reg.predict(X_test)))
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



##############################################
# Check overfit on the 'top feature XGB model'
##############################################



#####
# AUC plot for XGBoost
#####
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


#######
# AUC bar plots for slides
########

barWidth = 0.25
fig, ax = plt.subplots(figsize =(12, 8)) 
FM = [88, 86, 85, 79, 78] 
TF = [84, 84, 80, 72, 72] 
br1 = np.arange(len(FM)) 
br2 = [x + barWidth for x in br1] 
#F0E442
ax.bar(br1, FM, color ='#E69F00', width = barWidth, edgecolor ='grey', label ='Full 263 feature model') 
ax.bar(br2, TF, color ='#56B4F9', width = barWidth, edgecolor ='grey', label ='Top 25 feature model (based on RF perm.)') 
ax.set_xlabel('Machine learning method', fontweight ='bold', fontsize = 15) 
ax.set_ylabel('Performance (AUC)', fontweight ='bold', fontsize = 15) 
ax.set_xticks([r + barWidth/2 for r in range(len(FM))])
ax.set_xticklabels(['XGBoost', 'Random Forest', 'Multi-layer\nperceptron\n(untuned)', 'Support Vector\nMachine (untuned)', 'Logistic regression'])
ax.set_ylim(70, 90)
ax.set_yticks([70, 75, 80, 85, 90])
ax.legend()
plt.show()