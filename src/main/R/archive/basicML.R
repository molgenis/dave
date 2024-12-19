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

# Train model on ALL features
all_rf <-randomForest(ann_classificationVKGL~., data=all_train, ntree=1000, keep.forest=TRUE, importance=TRUE, do.trace=TRUE)

# Get most important features and generalize model
feat_imp_df <- importance(all_rf) %>% data.frame() %>% mutate(feature = row.names(.)) 
feat_imp_df[order(feat_imp_df$MeanDecreaseAccuracy, decreasing = T),]
# Plot: varImpPlot(fs1_rf)
topFeat <- feat_imp_df[order(feat_imp_df$MeanDecreaseAccuracy, decreasing = T)[1:25],]$feature
fs1_train <- all_train %>% select(c("ann_classificationVKGL", topFeat))
fs1_rf <-randomForest(ann_classificationVKGL~., data=fs1_train, ntree=1000, keep.forest=TRUE, importance=TRUE, do.trace=TRUE)

# Persist / load model for next time
#save(fs1_rf, file=paste(rootDir, "models", "rf_top25_deltaonly.Rdata", sep="/"))
load(paste(rootDir, "models", "rf_top25_deltaonly.Rdata", sep="/"))

# Apply model on test data
testData <- all_test # all_test, secr_test, memb_test, intr_test
rf.pr = predict(fs1_rf,type="prob",newdata=testData)[,2]
rf.pred = prediction(rf.pr, testData$ann_classificationVKGL)
rf.perf = performance(rf.pred, "tpr", "fpr")
auc <- performance(rf.pred,"auc")
auc <- unlist(slot(auc, "y.values"))
plot(rf.perf,main=paste0("Variant classification on peptide and\nfolding properties, RF ROC curve (AUC ",round(auc,2),")"),col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
auc # 97% total, 98% for secr and memb, 95% for intr

# PPV/NPV plot
rf.pv = performance(rf.pred, "ppv", "npv")
plot(rf.pv)
# Alternative
ppvnpv = performance(rf.pred, "ppv", "npv")
plot(ppvnpv,main="PPV/NPV plot",col=2,lwd=2)


###########################
# AlphaMissense performance only - using the same method as before
##################################
testData <- all_test # all_test, secr_test, memb_test, intr_test
am.pred <- prediction(testData$ann_am_pathogenicity, testData$ann_classificationVKGL)
am.perf = performance(am.pred, "tpr", "fpr")
auc <- performance(am.pred,"auc")
auc <- unlist(slot(auc, "y.values"))
auc

################
# Find good PPV/NPV thresholds using the test data draw on full data
################
data <- read.csv(freeze4)
rownames(data) <- paste0(data$gene, "/", data$UniProtID, ":", data$delta_aaSeq)
data <- subset(data, ann_classificationVKGL == "LP" | ann_classificationVKGL == "LB")
data <- data[!is.na(data$ann_am_pathogenicity), ]
data <- data[, colSums(data != 0) > 0]
data$ann_classificationVKGL <- as.factor(data$ann_classificationVKGL)
data$ann_proteinLocalization <- as.factor(data$ann_proteinLocalization)
full_test_data <- data[!all_draw, ]
test.pr = predict(fs1_rf,type="prob",newdata=full_test_data)[,2]
full_test_data <- cbind(full_test_data, test.pr)
ggplot(full_test_data, aes(test.pr, colour = ann_classificationVKGL, fill = ann_classificationVKGL)) +
  theme_classic() +
  geom_density(adjust=0.5, alpha = 0.25) +
  scale_fill_manual(name = "Classification", values = c("LB" = "green","LP" = "red")) +
  scale_colour_manual(name = "Classification", values = c("LB" = "green","LP" = "red"))
cutoff <- 0.75
tp <- sum(full_test_data[full_test_data$ann_classificationVKGL=="LP",'ann_am_pathogenicity'] >= cutoff)
fp <- sum(full_test_data[full_test_data$ann_classificationVKGL=="LB",'ann_am_pathogenicity'] >= cutoff)
tn <- sum(full_test_data[full_test_data$ann_classificationVKGL=="LB",'ann_am_pathogenicity'] < cutoff)
fn <- sum(full_test_data[full_test_data$ann_classificationVKGL=="LP",'ann_am_pathogenicity'] < cutoff)
ppv <-  100*tp/(tp+fp)
npv <-  100*tn/(tn+fn)
sens <- 100*tp/(tp+fn)
spec <- 100*tn/(fp+tn)
cat(paste("At cutoff ",round(cutoff,5)," we find PPV ",round(ppv),"%, NPV ",round(npv),"%, sensitivity ",round(sens),"% and specificity ",round(spec),"%.\n",sep=""))
# At cutoff 0.5 we find PPV 78%, NPV 94%, sensitivity 80% and specificity 93%.
# At cutoff 0.75 we find PPV 86%, NPV 90%, sensitivity 66% and specificity 97%.
# At cutoff 0.9 we find PPV 91%, NPV 88%, sensitivity 53% and specificity 98%.
# At cutoff 0.95 we find PPV 93%, NPV 85%, sensitivity 43% and specificity 99%.

#####
# Check for bias for mutant variable(s) that may correlate with genes with skewed LB/LP proportion
#####
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


#####
# VUS classification using model
#####
rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability"
freeze4 <- paste(rootDir, "data", "freeze4.csv.gz", sep="/")
vusPred <- read.csv(freeze4)
# For some reproduciblity, set a random seed, though RF is non-deterministic by design
set.seed(222)
# Nice row names
rownames(vusPred) <- paste0(vusPred$gene, "/", vusPred$UniProtID, ":", vusPred$delta_aaSeq)
# Select LB, LP and VUS, rows with AM, remove 0 cols
vusPred <- subset(vusPred, ann_classificationVKGL != "CF")
vusPred <- vusPred[!is.na(vusPred$ann_am_pathogenicity), ]
vusPred <- vusPred[, colSums(vusPred != 0) > 0]
# Apply model
load(paste(rootDir, "models", "rf_top25_deltaonly.Rdata", sep="/"))
predictions = predict(fs1_rf,type="prob",newdata=vusPred)[,2]
vusPred <- cbind(vusPred, predictions)
ggplot(vusPred, aes(predictions, colour = ann_classificationVKGL, fill = ann_classificationVKGL)) +
  theme_classic()+
  geom_density(adjust=0.5, alpha = 0.25)
vusOnlyPred <- subset(vusPred, ann_classificationVKGL == "VUS")
vusOnlyPredSecr <- subset(vusOnlyPred, ann_proteinLocalization == "membrane")
vusOnlyPredSecrLP <- subset(vusOnlyPredSecr, predictions >= 0.75)
nrow(vusOnlyPredSecrLP)
plot(vusOnlyPredSecrLP$ann_am_pathogenicity ~ vusOnlyPredSecrLP$predictions)
