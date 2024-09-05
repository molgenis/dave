library(xgboost)
library(caTools)
library(dplyr)
library(caret)

rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability"
freeze2 <- paste(rootDir, "data", "freeze2.csv.gz", sep="/")
data <- read.csv(freeze2)

####################
# Data preparation #
####################
# Nice row names
rownames(data) <- paste0(data$gene, "/", data$UniProtID, ":", data$delta_aaSeq)
# Select only LB and LP
data <- subset(data, ann_classificationVKGL == "LP" | ann_classificationVKGL == "LB")
# Remove all columns with only 0 values
data <- data[, colSums(data != 0) > 0]
# Select all columns with relevant factors or numerical variables for analysis
# We drop mutant and WT information here and only focus on the deltas
data <- data %>% select(contains(c("ann_classificationVKGL", "ann_proteinIschaperoned", "delta_", "mutant_")))
data <- data %>% select(-contains(c("ann_mutant_energy_SD", "_aaSeq")))
# Factorize categoricals
data$ann_classificationVKGL <- as.factor(data$ann_classificationVKGL)
# Check column types
sapply(data, class)

set.seed(42)
sample_split <- sample.split(Y = data$ann_classificationVKGL, SplitRatio = 0.7)
train_set <- subset(x = data, sample_split == TRUE)
test_set <- subset(x = data, sample_split == FALSE)

y_train <- as.integer(train_set$ann_classificationVKGL) - 1
y_test <- as.integer(test_set$ann_classificationVKGL) - 1
X_train <- train_set %>% select(-ann_classificationVKGL)
X_test <- test_set %>% select(-ann_classificationVKGL)

xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)
#xgb_params <- list(
#  booster = "gbtree",
#  eta = 0.01,
#  max_depth = 8,
#  gamma = 4,
#  subsample = 0.75,
#  colsample_bytree = 1,
#  objective = "multi:softprob",
#  eval_metric = "mlogloss",
#  num_class = length(levels(data$ann_classificationVKGL))
#)
xgb_params <- list(
  booster = "gbtree",
  eta = 0.05,  # Slightly higher learning rate for faster convergence
  max_depth = 6,  # Shallower trees to avoid overfitting
  gamma = 1,  # Moderate regularization
  subsample = 0.8,  # Slightly increased to leverage more data
  colsample_bytree = 0.8,  # Reduced to avoid overfitting
  objective = "binary:logistic",  # Corrected for binary classification
  eval_metric = "logloss"  # Changed metric to logloss
)

xgb_model <- xgb.train(
  params = xgb_params,
  data = xgb_train,
  nrounds = 5000,
  verbose = 1
)
xgb_model

xgb_preds <- predict(xgb_model, as.matrix(X_test), reshape = TRUE)
xgb_preds <- as.data.frame(xgb_preds)
colnames(xgb_preds) <- levels(data$ann_classificationVKGL)
xgb_preds

xgb_preds$PredictedClass <- apply(xgb_preds, 1, function(y) colnames(xgb_preds)[which.max(y)])
xgb_preds$ActualClass <- levels(data$ann_classificationVKGL)[y_test + 1]
xgb_preds

accuracy <- sum(xgb_preds$PredictedClass == xgb_preds$ActualClass) / nrow(xgb_preds)
accuracy

confusionMatrix(factor(xgb_preds$ActualClass), factor(xgb_preds$PredictedClass))

cv.res <- xgb.cv(data = xgb_train, nfold = 3, nrounds = 100, verbose = FALSE, objective = 'binary:logistic', eval_metric = 'auc', prediction = T)
it = which.max(cv.res$evaluation_log$test_auc_mean)
best.iter = cv.res$evaluation_log$iter[it]
plot(pROC::roc(response = y_train, predictor = cv.res$pred, levels=c(0, 1)), lwd=1.5) 
