library(randomForest)
library(caret)
library(pROC)
library(ROCR)
library(randomForestExplainer)

# interesting? https://jtr13.github.io/cc21fall2/introduction-to-xai-explainable-ai-in-r.html

rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability"
freeze1 <- paste(rootDir, "data", "freeze1.csv", sep="/")
data <- read.csv(freeze1)
data <- subset(data, ann_classificationVKGL == "LP" | ann_classificationVKGL == "LB")
data$ann_classificationVKGL <- as.factor(data$ann_classificationVKGL)

set.seed(222)
draw <- sample(c(TRUE, FALSE), nrow(data), replace=TRUE, prob=c(0.8, 0.2))
train <- data[draw, ]
test <- data[!draw, ]

# worked before not not now? optimal mtry for RF
#bestmtry <- tuneRF(train,train$classificationVKGL,stepFactor = 1.2, improve = 0.01, trace=T, plot= T) 

rf <-randomForest(ann_classificationVKGL~., data=train, mtry=7, ntree=1000, keep.forest=TRUE, importance=TRUE, xtest=subset(test, select=-ann_classificationVKGL))
rf.pr = predict(rf,type="prob",newdata=test)[,2]
rf.pred = prediction(rf.pr, test$ann_classificationVKGL)
rf.perf = performance(rf.pred, "tpr", "fpr")
auc <- performance(rf.pred,"auc")
auc <- unlist(slot(auc, "y.values"))
plot(rf.perf,main=paste0("Variant classification on basic protein properties, RF ROC curve (AUC ",round(auc,2),")"),col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
auc

rf.perf.pnpv = performance(rf.pred, "ppv", "npv")
plot(rf.perf.pnpv,main=paste0("Variant classification on basic protein properties, RF ROC curve (AUC ",round(auc,2),")"),col=2,lwd=2)

# from https://datasciencechronicle.wordpress.com/2014/03/17/r-tips-part2-rocr-example-with-randomforest/

# feature importance?
varImpPlot(rf)
measure_importance(rf)

# PCA?
# affinity-prop clustering?
