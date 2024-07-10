library(randomForest)
library(caret)
library(pROC)
library(ROCR)

rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability"
freeze1 <- paste(rootDir, "data", "freeze1.csv", sep="/")
data <- read.csv(freeze1)
data <- subset(data, classificationVKGL == "LP" | classificationVKGL == "LB")
data$classificationVKGL <- as.factor(data$classificationVKGL)

set.seed(222)
draw <- sample(c(TRUE, FALSE), nrow(data), replace=TRUE, prob=c(0.8, 0.2))
train <- data[draw, ]
test <- data[!draw, ]

bestmtry <- tuneRF(train,train$classificationVKGL,stepFactor = 1.2, improve = 0.01, trace=T, plot= T) 

rf <-randomForest(classificationVKGL~., data=train, mtry=7, ntree=1000, keep.forest=TRUE, importance=TRUE, xtest=subset(test, select=-classificationVKGL))
rf.pred = prediction(rf.pr, test$classificationVKGL)
rf.perf = performance(rf.pred,"tpr","fpr")
auc <- performance(rf.pred,"auc")
auc <-unlist(slot(auc, "y.values"))
auc <- round(auc,2)
plot(rf.perf,main=paste0("Variant classification on basic protein properties, RF ROC curve (AUC ",auc,")"),col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")



# PCA?
# affinity-prop clustering?
