library(randomForest)
library(caret)
library(pROC)
library(ROCR)
library(randomForestExplainer)
library(dplyr)

# interesting? https://jtr13.github.io/cc21fall2/introduction-to-xai-explainable-ai-in-r.html

rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability"
freeze4 <- paste(rootDir, "data", "freeze4.csv.gz", sep="/")
data <- read.csv(freeze4)

####################
# Data preparation #
####################
# For some reproduciblity, set a random seed, though RF is non-deterministic by design
set.seed(222)
# Nice row names
rownames(data) <- paste0(data$gene, "/", data$UniProtID, ":", data$delta_aaSeq)
# Select only LB and LP
data <- subset(data, ann_classificationVKGL == "LP" | ann_classificationVKGL == "LB")
# For AM result, remove rows for which ann_am_pathogenicity = NA
data <- data[!is.na(data$ann_am_pathogenicity), ]
# Remove all columns with only 0 values
data <- data[, colSums(data != 0) > 0]
# Select all columns with relevant factors or numerical variables for analysis
# We drop mutant and WT information here and only focus on the deltas
# With or without "ann_am_pathogenicity"
data <- data %>% select(contains(c("ann_classificationVKGL", "ann_proteinIschaperoned", "ann_am_pathogenicity", "ann_proteinLocalization", "delta_")))
data <- data %>% select(-contains(c("ann_mutant_energy_SD", "_aaSeq")))

# Factorize categoricals
data$ann_classificationVKGL <- as.factor(data$ann_classificationVKGL)
data$ann_proteinLocalization <- as.factor(data$ann_proteinLocalization)

# Check column types
sapply(data, class)

# Divide and subset
propTrain <- 0.8
propTest <- 0.2

all_draw <- sample(c(TRUE, FALSE), nrow(data), replace=TRUE, prob=c(propTrain, propTest))
all_train <- data[all_draw, ]
all_test <- data[!all_draw, ]

secr_test <- subset(all_test, ann_proteinLocalization == "secreted")
memb_test <- subset(all_test, ann_proteinLocalization == "membrane")
intr_test <- subset(all_test, ann_proteinLocalization == "intracellular")

# Train complete model
all_rf <-randomForest(ann_classificationVKGL~., data=all_train, ntree=1000, keep.forest=TRUE, importance=TRUE, do.trace=TRUE)

# Get most important features and generalize model
feat_imp_df <- importance(all_rf) %>% data.frame() %>% mutate(feature = row.names(.)) 
feat_imp_df[order(feat_imp_df$MeanDecreaseAccuracy, decreasing = T),]
topFeat <- feat_imp_df[order(feat_imp_df$MeanDecreaseAccuracy, decreasing = T)[1:25],]$feature
fs1_train <- all_train %>% select(c("ann_classificationVKGL", topFeat))
fs1_rf <-randomForest(ann_classificationVKGL~., data=fs1_train, ntree=1000, keep.forest=TRUE, importance=TRUE, do.trace=TRUE)
#save(fs1_rf, file="rf_top25_deltaonly.Rdata")
#load(paste(rootDir, "data", "rf_top25_deltaonly.Rdata", sep="/"))
rf <- fs1_rf # all_rf, secr_rf, memb_intr_rf
testData <- all_test # all_test, secr_test, memb_test, intr_test
rf.pr = predict(rf,type="prob",newdata=testData)[,2]
rf.pred = prediction(rf.pr, testData$ann_classificationVKGL)
rf.perf = performance(rf.pred, "tpr", "fpr")
auc <- performance(rf.pred,"auc")
auc <- unlist(slot(auc, "y.values"))
plot(rf.perf,main=paste0("Variant classification on peptide and\nfolding properties, RF ROC curve (AUC ",round(auc,2),")"),col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
auc # 97% total, 98% for secr and memb, 95% for intr

rf.pv = performance(rf.pred, "ppv", "npv")
plot(rf.pv)


# Extra stuff: PPV/NPV, feature importance, etc
ppvnpv = performance(rf.pred, "ppv", "npv")
plot(ppvnpv,main="PPV/NPV plot",col=2,lwd=2)
# from https://datasciencechronicle.wordpress.com/2014/03/17/r-tips-part2-rocr-example-with-randomforest/
# feature importance?
varImpPlot(rf)
measure_importance(rf)


# Check for bias for mutant variable(s) that may correlate with genes with skewed LB/LP proportion
biasCheck <- read.csv(freeze4)
biasCheckGeneRes <- data.frame()
uniqGenes <- unique(biasCheck$gene)
rownames(biasCheck) <- paste0(biasCheck$gene, "/", biasCheck$UniProtID, ":", biasCheck$delta_aaSeq)
biasCheck <- subset(biasCheck, ann_classificationVKGL == "LP" | ann_classificationVKGL == "LB")
for(selectGene in uniqGenes){
  cat(paste0(selectGene,"\n"))
  byGene <- subset(biasCheck, gene == selectGene)
  lpS <- sum(str_count(byGene$ann_classificationVKGL, "LP"))
  lbS <- sum(str_count(byGene$ann_classificationVKGL, "LB"))
  if(lpS == 0 & lbS > 0){
    propLP <- 0
  }else if(lbS == 0 & lpS > 0){
    propLP <- 1 
  }else{
    propLP <- lpS/(lpS+lbS)
  }
  geneRes <- data.frame(propLP = propLP, WT_Xc2.lambda.28 = byGene$WT_Xc2.lambda.28[1])
  biasCheckGeneRes <- rbind(biasCheckGeneRes, geneRes)
}
plot(biasCheckGeneRes$propLP ~ biasCheckGeneRes$WT_Xc2.lambda.28)
cor.test(biasCheckGeneRes$propLP, biasCheckGeneRes$WT_Xc2.lambda.28)
cor.test(biasCheckGeneRes$propLP, biasCheckGeneRes$WT_Xc2.lambda.28, method="kendall")
# Indeed: genes with a low proportion of LP also have low WT_Xc2.lambda.28
# Therefore we have detected a bias


