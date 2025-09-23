library(crunch)     # To compress results
library(R.utils)    # for 'gunzip', 'mkdirs'
library(dplyr)      # Data manipulation
library(tidyr)      # Data manipulation
library(xgboost)
library(caret)
library(randomForest)
library(reshape2)
library(pheatmap)
library(pROC)
library(rstanarm)
library(corrplot)
library(fastshap)
library(ggplot2)


#######################
# Adjustable settings #
#######################
rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability" # root directory that contains README.md, data/, img/, out/, src/, etc.


###################################
# Invariants: seperators and seed #
###################################
CONCAT_SEP <- "~" # used to combine UniProt hierarchy: category -> type -> description -> etc
VALUE_SEP <- "|" # used to combine multiple annotation values in the output data freeze
set.seed(1234)


##################################
# Load data and helper functions #
##################################
frz6loc <- paste(rootDir, "data", "freeze6.csv.gz", sep="/")
frz6 <- read.csv(frz6loc)
source(paste(rootDir, "src", "main", "R", "SHAPDecisionPlot.R", sep="/"))


############################################
# Define which features will be target, 'explainers', or removed.
#######################################


predictionTarget <- "ann_classificationVKGL"

modelFeatures <- c(
  # from GeoNet
  "delta_DNAs_cumu_bin", "delta_RNAs_cumu_bin", "delta_ProtS_cumu_bin",
  # from GLM-Score
  "delta_DNAb_binding_affinity_pKd", "delta_RNAb_binding_affinity_pKd",
  # from P2Rank
  "delta_ligand_nr_of_predicted_pockets", "delta_ligand_rank1_sas_points",
  # from R-Peptides
  "delta_charge", "delta_hydrophobicMoment", "delta_hydrophobicity", "delta_isoElecPoint",
  # from FoldX
  "delta_total.energy"
)

variantContext <- c("gene",
                "TranscriptID",
                "UniProtID",
                "dna_variant_chrom",
                "dna_variant_pos",
                "dna_variant_ref",
                "dna_variant_alt",
                "dna_variant_assembly",
                "delta_aaSeq",
                "ann_proteinIschaperoned",
                "ann_proteinLocalization",
                "ann_am_pathogenicity",
                "seqFt",
                "chain")


##########################################################
# Clean, slice, investigate features and train submodels #
##########################################################

# remove features that are never used
keep <- c(predictionTarget, modelFeatures, variantContext)
frz6 <- frz6 %>% select(any_of(keep))
# divide into train/test and new observations
frz6_train_test <- frz6 %>% filter(ann_classificationVKGL %in% c("LB", "LP"))
frz6_new_observ <- frz6 %>% filter(ann_classificationVKGL %in% c("VUS")) # no 'CF'
# factorize now, to get LB/LP and VUS/CF factors per set
frz6_train_test$ann_classificationVKGL <- as.factor(frz6_train_test$ann_classificationVKGL)
frz6_new_observ$ann_classificationVKGL <- as.factor(frz6_new_observ$ann_classificationVKGL)
# remove context from train/test, divide into 80% train / 20% test
frz6_train_test <- frz6_train_test %>% select(-all_of(variantContext))
train_idx <- createDataPartition(frz6_train_test$ann_classificationVKGL, p = 0.8, list = FALSE)
frz6_train <- frz6_train_test[train_idx,]
frz6_test <- frz6_train_test[-train_idx,]

# calculate R-square correlation across features and plot
cor_r2_mat <- cor(frz6_train[ , !(names(frz6_train) %in% predictionTarget) ])^2
corrplot(cor_r2_mat,  type = "upper", tl.cex=0.7, tl.col = "black", method = "color", addCoef.col = "black",  number.cex = 0.5)

# fit basic Random Forest model, show feature importance and AUC
rf_model <- randomForest(ann_classificationVKGL ~ .,data = frz6_train, importance = TRUE)
preds <- predict(rf_model, newdata = frz6_test)
print(rf_model)
importance(rf_model) # optional plot: varImpPlot(rf_model)
rf_probs <- predict(rf_model, newdata = frz6_test, type = "prob")
roc_obj <- roc(frz6_test$ann_classificationVKGL, rf_probs[,"LP"])
auc(roc_obj) # optional plot: plot(roc_obj, col = "blue", lwd = 2, main = "Random Forest ROC Curve")



####
# Apply RF model on new observations to get predictions
####

frz6_new_observ_pred_prob <- predict(rf_model, newdata = frz6_new_observ, type = "prob")  # class probabilities
#aa <- cbind(frz6_new_observ, rf_pred_prob)

# Wrap the predict function (needs to return numeric probs for SHAP)
pfun <- function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, "LP"]  # pick probability for "LP"
}
frz6_new_observ_only_modelFeatures<- frz6_new_observ %>% select(any_of(modelFeatures))
# nsim should always be >1 to obtain baseline,  but more is better
frz6_train_without_target <- frz6_train %>% select(-all_of(predictionTarget))
frz6_new_observ_explain <- fastshap::explain(rf_model, shap_only = FALSE,
  X = frz6_train_without_target, newdata = frz6_new_observ_only_modelFeatures, pred_wrapper = pfun, nsim = 10, adjust = TRUE
)
# 
frz6_new_observ_shap <- as.data.frame(frz6_new_observ_explain$shapley_values)
colnames(frz6_new_observ_shap) <- paste0(colnames(frz6_new_observ_explain$shapley_values), ".sph")
# Add columns for base and final probabilities
frz6_new_observ_shap$BaseProbability.sph  <- frz6_new_observ_explain$baseline
frz6_new_observ_shap$FinalProbability.sph <- rowSums(frz6_new_observ_shap) # + frz6_new_observ_explain$baseline
plot(sort(frz6_new_observ_shap$FinalProbability.sph))

# merge back with context for plot!
all_vus <- cbind(frz6_new_observ, frz6_new_observ_shap)
all_vus <- cbind(all_vus, frz6_new_observ_pred_prob)
#Sanity check
plot(all_vus$FinalProbability.sph, all_vus$LP)

all_vus_sorted <- all_vus %>% arrange(FinalProbability.sph)

write.csv(all_vus_sorted[c(1:10, (nrow(all_vus_sorted)-9):(nrow(all_vus_sorted))),], file="vus_top10_bottom10.csv")

shapDecisionPlot(all_vus_sorted[1,])
shapDecisionPlot(all_vus_sorted[2,])
shapDecisionPlot(all_vus_sorted[3,])

shapDecisionPlot(all_vus_sorted[11367,])
shapDecisionPlot(all_vus_sorted[11368,])
shapDecisionPlot(all_vus_sorted[11369,])

##
# Explain predictions using SHAP heuristic
##


#####
# Combine all results and save to final predictions and explanations
####

