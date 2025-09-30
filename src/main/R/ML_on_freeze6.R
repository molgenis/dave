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
library(cutpointr)
library(vcfR)


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
vus_pred_loc <- paste(rootDir, "data", "vkgl_apr2024_VUS_pred.csv.gz", sep="/")


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
corrFeatRelabel <- read.csv(paste(rootDir, "data", "12features.csv", sep="/"))
corrFeatRelabel$Name <- sub("\\*+$", "", corrFeatRelabel$Name)
#feat$Name <- ifelse(nzchar(feat$Unit), paste0(feat$Name, " (in ", feat$Unit, ")"), feat$Name) # adds the unit, but not needed
old_names <- rownames(cor_r2_mat) # Get the current names from cor_r2_mat
lookup <- setNames(corrFeatRelabel$Name, corrFeatRelabel$Feature) # Create a lookup vector from feat
new_names <- lookup[old_names] # Replace row and column names by matching
rownames(cor_r2_mat) <- new_names
colnames(cor_r2_mat) <- new_names
pdf_plot_loc <- paste(rootDir, "img", "corrplot.pdf", sep="/")
cairo_pdf(file = pdf_plot_loc, width = 5, height = 5)  # adjust size as needed
corrplot(cor_r2_mat,  type = "upper", tl.cex=0.7, tl.col = "black", method = "color", addCoef.col = "black",  number.cex = 0.5)
dev.off()

# fit basic Random Forest model, show feature importance and AUC
rf_model <- randomForest(ann_classificationVKGL ~ .,data = frz6_train, importance = TRUE)
preds <- predict(rf_model, newdata = frz6_test)
print(rf_model)
importance(rf_model) # optional plot: varImpPlot(rf_model)
rf_probs <- predict(rf_model, newdata = frz6_test, type = "prob")
roc_obj <- pROC::roc(frz6_test$ann_classificationVKGL, rf_probs[,"LP"])
pROC::auc(roc_obj) # optional plot: plot(roc_obj, col = "blue", lwd = 2, main = "Random Forest ROC Curve")



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
frz6_new_observ_shap$FinalProbability.sph <- rowSums(frz6_new_observ_shap)
# Show predicted probabilities
plot(sort(frz6_new_observ_shap$FinalProbability.sph))

# Merge back with context for further plots
all_vus <- cbind(frz6_new_observ, frz6_new_observ_shap)
all_vus <- cbind(all_vus, frz6_new_observ_pred_prob)
# Sanity check: are model probs and SHAP probs lined up correctly?
plot(all_vus$FinalProbability.sph, all_vus$LP, xlim=c(0,1), ylim=c(0,1))
# Sanity check: plot vs AlphaMissense (notice imputed value around 0.39)
plot(all_vus$ann_am_pathogenicity, all_vus$LP, xlim=c(0,1), ylim=c(0,1))
# Sort by final prob
all_vus_sorted <- all_vus %>% arrange(FinalProbability.sph)
# Write to disk for later fast reuse
write.csv.gz(all_vus_sorted, file=vus_pred_loc, row.names=F)

##########
# Start point for VUS analysis after writing results previously
##########

# Load VUS predictions back in
all_vus_sorted <- read.csv(file=vus_pred_loc)
# Plot 'most benign' and 'most pathogenic' predictions
plotrows <- c(1,2,3,4,5,11217,11218,11219,11220,11221)
for(plotrow in plotrows)
{
  row <- all_vus_sorted[plotrow,]
  p <- shapDecisionPlot(row)
  pdf_plot_loc <- paste(rootDir, "img", paste0(row$gene, "_", row$delta_aaSeq, ".pdf"), sep="/")
  png_plot_loc <- paste(rootDir, "img", paste0(row$gene, "_", row$delta_aaSeq, ".png"), sep="/")
  ggsave(filename = pdf_plot_loc, plot = p, device = cairo_pdf, width = 10, height = 4) # height was 6.25
  #ggsave(filename = png_plot_loc, plot = p, device = "png", width = 9.777778, height = 5.5) # 16:9 as PNG for full screen slides
}
# as table
all_vus_sorted$dna <- paste0(all_vus_sorted$dna_variant_chrom,":",all_vus_sorted$dna_variant_pos," ",all_vus_sorted$dna_variant_ref,">",all_vus_sorted$dna_variant_alt)
all_vus_sorted[plotrows,c("gene","UniProtID","dna","delta_aaSeq","LP")]
# find with affected ligand top pocket
#all_vus_ligand_aff <- all_vus_sorted %>% arrange(delta_ligand_rank1_sas_points)
#all_vus_ligand_aff[plotrows,c("gene","UniProtID","dna","delta_aaSeq","LP","delta_ligand_rank1_sas_points")]
#p <- shapDecisionPlot(all_vus_ligand_aff[11221,])

# Merge with variants that can received a classification in the meantime
fr6_from_vus_to_lp_lb_loc <- paste(rootDir, "data", "fr6-vkgl-clf-vus-to-lp-lb-apr2024-july2025.csv", sep="/")
fr6_from_vus_to_lp_lb <- read.csv(fr6_from_vus_to_lp_lb_loc)
vus_changed <- merge(x = fr6_from_vus_to_lp_lb, y = all_vus_sorted, by.x = c(variantContext,modelFeatures,"ann_classificationVKGL"), by.y = c(variantContext,modelFeatures,"ann_classificationVKGL"))
# New classification vs prediction
plot(as.factor(vus_changed$new_classification), vus_changed$FinalProbability.sph)
# Make SHAP breakdown plot of each prediction
for(i in 1:nrow(vus_changed)) {
  row <- vus_changed[i,]
  p <- shapDecisionPlot(row)
  pdf_plot_loc <- paste(rootDir, "img", paste0(row$gene, "_", row$delta_aaSeq, ".pdf"), sep="/")
  ggsave(filename = pdf_plot_loc, plot = p, device = cairo_pdf, width = 10, height = 4)
}
# as table
vus_changed_sorted <- vus_changed %>% arrange(LP)
vus_changed_sorted$verdict <- ifelse(vus_changed_sorted$LP >= 0.286,"P","B")# based on 0.286, see below
vus_changed_sorted$dna <- paste0(vus_changed_sorted$dna_variant_chrom,":",vus_changed_sorted$dna_variant_pos," ",vus_changed_sorted$dna_variant_ref,">",vus_changed_sorted$dna_variant_alt)
vus_changed_sorted[,c("gene","TranscriptID","UniProtID","dna","delta_aaSeq","LP", "verdict","new_classification")]


#### Now on ClinVar data
# download from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20250923.vcf.gz
# or later from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2025/clinvar_20250923.vcf.gz
clinvar_loc <- paste(rootDir, "data", "clinvar_20250923.vcf.gz", sep="/")
clinvarVCF <- read.vcfR(clinvar_loc)
clinvar <- as.data.frame(clinvarVCF@fix)
vus_changed_clinv <- merge(x = clinvar, y = all_vus_sorted, by.x = c("CHROM", "POS", "REF", "ALT"), by.y = c( "dna_variant_chrom", "dna_variant_pos", "dna_variant_ref", "dna_variant_alt"))
vus_changed_clinv_LP <- subset(vus_changed_clinv, grepl("CLNSIG=(Likely_pathogenic|Pathogenic)", INFO))
vus_changed_clinv_LB <- subset(vus_changed_clinv, grepl("CLNSIG=(Likely_benign|Benign)", INFO))
vus_changed_clinv_LP$new_classification <- "LP/P"
vus_changed_clinv_LB$new_classification <- "LB/B"
vus_changed_clinv_both <- rbind(vus_changed_clinv_LP, vus_changed_clinv_LB)
plot(as.factor(vus_changed_clinv_both$new_classification), vus_changed_clinv_both$FinalProbability.sph)
# find with affected ligand top pocket
vus_changed_clinv_both_ligand_aff <- vus_changed_clinv_both %>% arrange(delta_ligand_rank1_sas_points)
vus_changed_clinv_both_ligand_aff[c(1,2,3,707,708,709),c("gene","UniProtID","dna","delta_aaSeq","LP","delta_ligand_rank1_sas_points")]
p <- shapDecisionPlot(all_vus_ligand_aff[11221,])
p

# Determine optimal threshold on ClinVar using Youden's Index
cutpointDF <- subset(vus_changed_clinv_both, new_classification == "LB/B" | new_classification == "LP/P")
opt_cut <- cutpointr(cutpointDF, FinalProbability.sph, new_classification, direction = ">=", pos_class = "LP/P", neg_class = "LB/B", method = maximize_metric, metric = youden)
youdenIndex <- opt_cut$optimal_cutpoint # = 0.286, use as default threshold
tp <- sum(cutpointDF[cutpointDF$new_classification=="LP/P",'FinalProbability.sph'] >= youdenIndex)
fp <- sum(cutpointDF[cutpointDF$new_classification=="LB/B",'FinalProbability.sph'] >= youdenIndex)
tn <- sum(cutpointDF[cutpointDF$new_classification=="LB/B",'FinalProbability.sph'] < youdenIndex)
fn <- sum(cutpointDF[cutpointDF$new_classification=="LP/P",'FinalProbability.sph'] < youdenIndex)
ppv <- 100 *tp/(tp+fp)
npv <- 100 *tn/(tn+fn)
sens <- tp / (tp + fn)*100
spec <- tn / (tn + fp)*100
cat(paste("in ClinVar data we find", tp, "TP,", fp, "FP,", tn, "TN and", fn, "FN"))
cat(paste("this means ", ppv, "PPV,", npv, "NPV,", sens, "sens and", spec, "spec"))

# Apply this threshold on VKGL
#cutpointDF <- subset(vus_changed, new_classification == "LB" | new_classification == "LP")
#opt_cut <- cutpointr(cutpointDF, FinalProbability.sph, new_classification, direction = ">=", pos_class = "LP", neg_class = "LB", method = maximize_metric, metric = youden)
#youdenIndex <- opt_cut$optimal_cutpoint # here, 0.278, but we're using ClinVar's youden
youdenIndex <- 0.286
tp <- sum(vus_changed[vus_changed$new_classification=="LP",'FinalProbability.sph'] >= youdenIndex)
fp <- sum(vus_changed[vus_changed$new_classification=="LB",'FinalProbability.sph'] >= youdenIndex)
tn <- sum(vus_changed[vus_changed$new_classification=="LB",'FinalProbability.sph'] < youdenIndex)
fn <- sum(vus_changed[vus_changed$new_classification=="LP",'FinalProbability.sph'] < youdenIndex)
ppv <- 100 *tp/(tp+fp)
npv <- 100 *tn/(tn+fn)
sens <- tp / (tp + fn)*100
spec <- tn / (tn + fp)*100
cat(paste("when applied to VKGL we find", tp, "TP,", fp, "FP,", tn, "TN and", fn, "FN"))
cat(paste("this means ", ppv, "PPV,", npv, "NPV,", sens, "sens and", spec, "spec"))
# --> 4 TP, 2 FP, 5 TN and 1 FN
