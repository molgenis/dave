library(xgboost)
library(caTools)
library(dplyr)
library(caret)

rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability"
freeze1 <- paste(rootDir, "data", "freeze1.csv", sep="/")
data <- read.csv(freeze1)
data <- subset(data, classificationVKGL == "LP" | classificationVKGL == "LB")
data$classificationVKGL <- as.factor(data$classificationVKGL)
data$Pdb <- NULL
data$mutatedAAseq <- NULL
data$assembly <- NULL
data$chrom <- NULL
data$ref <- NULL
data$alt <- NULL
data$gene <- NULL
data$protChange <- NULL
data$transcript <- NULL
data$uniprot <- NULL
data$protType <- NULL

set.seed(42)
sample_split <- sample.split(Y = data$classificationVKGL, SplitRatio = 0.7)
train_set <- subset(x = data, sample_split == TRUE)
test_set <- subset(x = data, sample_split == FALSE)

y_train <- as.integer(train_set$classificationVKGL) - 1
y_test <- as.integer(test_set$classificationVKGL) - 1
X_train <- train_set %>% select(-classificationVKGL)
X_test <- test_set %>% select(-classificationVKGL)

xgb_train <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
xgb_test <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)
xgb_params <- list(
  booster = "gbtree",
  eta = 0.01,
  max_depth = 8,
  gamma = 4,
  subsample = 0.75,
  colsample_bytree = 1,
  objective = "multi:softprob",
  eval_metric = "mlogloss",
  num_class = length(levels(data$classificationVKGL))
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
colnames(xgb_preds) <- levels(data$classificationVKGL)
xgb_preds

xgb_preds$PredictedClass <- apply(xgb_preds, 1, function(y) colnames(xgb_preds)[which.max(y)])
xgb_preds$ActualClass <- levels(data$classificationVKGL)[y_test + 1]
xgb_preds

accuracy <- sum(xgb_preds$PredictedClass == xgb_preds$ActualClass) / nrow(xgb_preds)
accuracy

confusionMatrix(factor(xgb_preds$ActualClass), factor(xgb_preds$PredictedClass))

cv.res <- xgb.cv(data = xgb_train, nfold = 3, nrounds = 100, verbose = FALSE, objective = 'binary:logistic', eval_metric = 'auc', prediction = T)
it = which.max(cv.res$evaluation_log$test_auc_mean)
best.iter = cv.res$evaluation_log$iter[it]
plot(pROC::roc(response = y_train, predictor = cv.res$pred, levels=c(0, 1)), lwd=1.5) 
