library(randomForest)
library(caret)
library(pROC)
library(ROCR)
library(randomForestExplainer)
library(dplyr)

# interesting? https://jtr13.github.io/cc21fall2/introduction-to-xai-explainable-ai-in-r.html

rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability"
freeze2 <- paste(rootDir, "data", "freeze2.csv.gz", sep="/")
data <- read.csv(freeze2)

####################
# Data preparation #
####################
# For some reproduciblity, set a random seed, though RF is non-deterministic by design
set.seed(222)
# Nice row names
rownames(data) <- paste0(data$gene, "/", data$UniProtID, ":", data$delta_aaSeq)
# Select only LB and LP
data <- subset(data, ann_classificationVKGL == "LP" | ann_classificationVKGL == "LB")
# Remove all columns with only 0 values
data <- data[, colSums(data != 0) > 0]
# Select all columns with relevant factors or numerical variables for analysis
# We drop mutant and WT information here and only focus on the deltas
data <- data %>% select(contains(c("ann_classificationVKGL", "ann_proteinIschaperoned", "ann_proteinLocalization", "delta_", "mutant_")))
data <- data %>% select(-contains(c("ann_mutant_energy_SD", "_aaSeq")))

# Factorize categoricals
data$ann_classificationVKGL <- as.factor(data$ann_classificationVKGL)
data$ann_proteinLocalization <- as.factor(data$ann_proteinLocalization)

# Check column types
sapply(data, class)

# Subselect on localization, chaperonization, etc
#secr_data <- subset(data, ann_proteinLocalization == "secreted")
#memb_intr_data <- subset(data, ann_proteinLocalization == "membrane" | ann_proteinLocalization == "intracellular")
#data <- subset(data, ann_proteinIschaperoned == FALSE)

propTrain <- 0.8
propTest <- 0.2

all_draw <- sample(c(TRUE, FALSE), nrow(data), replace=TRUE, prob=c(propTrain, propTest))
all_train <- data[all_draw, ]
all_test <- data[!all_draw, ]

secr_train_data <- subset(all_train, ann_proteinLocalization == "secreted")
secr_train_draw <- sample(c(TRUE, FALSE), nrow(secr_train_data), replace=TRUE, prob=c(propTrain, propTest))
secr_train <- secr_train_data[secr_train_draw, ]

secr_test_data <- subset(all_test, ann_proteinLocalization == "secreted")
secr_test_draw <- sample(c(TRUE, FALSE), nrow(secr_test_data), replace=TRUE, prob=c(propTrain, propTest))
secr_test <- secr_test_data[!secr_test_draw, ]

memb_intr_train_data <- subset(all_train, ann_proteinLocalization == "membrane" | ann_proteinLocalization == "intracellular")
memb_intr_train_draw <- sample(c(TRUE, FALSE), nrow(memb_intr_train_data), replace=TRUE, prob=c(propTrain, propTest))
memb_intr_train <- memb_intr_train_data[memb_intr_train_draw, ]

memb_intr_test_data <- subset(all_test, ann_proteinLocalization == "membrane" | ann_proteinLocalization == "intracellular")
memb_intr_test_draw <- sample(c(TRUE, FALSE), nrow(memb_intr_test_data), replace=TRUE, prob=c(propTrain, propTest))
memb_intr_test <- memb_intr_test_data[!memb_intr_test_draw, ]

# worked before not not now? optimal mtry for RF
#bestmtry <- tuneRF(train,train$classificationVKGL,stepFactor = 1.2, improve = 0.01, trace=T, plot= T) # mtry=16,  mtry=round(sqrt(length(data))) ?
all_rf <-randomForest(ann_classificationVKGL~., data=all_train, ntree=100, keep.forest=TRUE, importance=TRUE, do.trace=TRUE)
secr_rf <-randomForest(ann_classificationVKGL~., data=secr_train, ntree=1000, keep.forest=TRUE, importance=TRUE, do.trace=TRUE)
memb_intr_rf <-randomForest(ann_classificationVKGL~., data=memb_intr_train, ntree=100, keep.forest=TRUE, importance=TRUE, do.trace=TRUE)

rf <- secr_rf # all_rf, secr_rf, memb_intr_rf
testData <- secr_test # all_test, secr_test, memb_intr_test
rf.pr = predict(rf,type="prob",newdata=testData)[,2]
rf.pred = prediction(rf.pr, testData$ann_classificationVKGL)
rf.perf = performance(rf.pred, "tpr", "fpr")
auc <- performance(rf.pred,"auc")
auc <- unlist(slot(auc, "y.values"))
plot(rf.perf,main=paste0("Variant classification on peptide and\nfolding properties, RF ROC curve (AUC ",round(auc,2),")"),col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
auc

rf.pv = performance(rf.pred, "ppv", "npv")
plot(rf.pv)


# Extra stuff: PPV/NPV, feature importance, etc
ppvnpv = performance(rf.pred, "ppv", "npv")
plot(ppvnpv,main="PPV/NPV plot",col=2,lwd=2)
# from https://datasciencechronicle.wordpress.com/2014/03/17/r-tips-part2-rocr-example-with-randomforest/
# feature importance?
varImpPlot(rf)
measure_importance(rf)
# PCA?
# affinity-prop clustering?
