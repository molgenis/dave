library(crunch)
library(dplyr)
library(caret)
library(randomForest)
library(pROC)
library(corrplot)
library(cutpointr)
library(vcfR)


#######################
# Adjustable settings #
#######################
rootDir <- "/Users/joeri/git/dave1" # root directory that contains README.md, data/, img/, out/, src/, etc.


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
pdf_plot_loc <- paste(rootDir, "img", "model-feature-correlations.pdf", sep="/")
cairo_pdf(file = pdf_plot_loc, width = 5, height = 5)  # adjust size as needed
corrplot(cor_r2_mat,  type = "upper", tl.cex=0.7, tl.col = "black", method = "color", addCoef.col = "black",  number.cex = 0.5)
dev.off()

# fit basic Random Forest model, show feature importance and AUC
rf_model <- randomForest(ann_classificationVKGL ~ .,data = frz6_train, importance = TRUE)
frz6_test_pred_prob <- predict(rf_model, newdata = frz6_test, type = "prob")  # class probabilities
# determine optimal prob threshold based on test data
comb_test_pred <- cbind(frz6_test, frz6_test_pred_prob)
plot(comb_test_pred$ann_classificationVKGL, comb_test_pred$LP)
opt_cut <- cutpointr(comb_test_pred, LP, ann_classificationVKGL, direction = ">=", pos_class = "LP", neg_class = "LB", method = maximize_metric, metric = youden)
youdenIndex <- opt_cut$optimal_cutpoint
youdenIndex # optimal threshold is 0.286
# print model info and AUC
print(rf_model)
importance(rf_model) # optional plot: varImpPlot(rf_model)
roc_obj <- pROC::roc(comb_test_pred$ann_classificationVKGL, comb_test_pred$LP)
pROC::auc(roc_obj) # optional plot: plot(roc_obj, col = "blue", lwd = 2, main = "Random Forest ROC Curve")
#save(rf_model,file = paste(rootDir, "models", "dave1_rf_model.RData", sep="/"))



################
# Apply RF model on new observations to get predictions and use SHAP to explain
################
# New predictions
frz6_new_observ_pred_prob <- predict(rf_model, newdata = frz6_new_observ, type = "prob")  # class probabilities
# To get SHAP, first wrap the predict function (needs to return numeric probs for SHAP)
pfun <- function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, "LP"]
}
frz6_new_observ_only_modelFeatures<- frz6_new_observ %>% select(any_of(modelFeatures))
frz6_train_without_target <- frz6_train %>% select(-all_of(predictionTarget))
frz6_new_observ_explain <- fastshap::explain(rf_model, shap_only = FALSE,
  X = frz6_train_without_target, newdata = frz6_new_observ_only_modelFeatures, pred_wrapper = pfun, nsim = 10, adjust = TRUE
)
frz6_new_observ_shap <- as.data.frame(frz6_new_observ_explain$shapley_values)
colnames(frz6_new_observ_shap) <- paste0(colnames(frz6_new_observ_explain$shapley_values), ".sph")
# Add columns for base and final probabilities
frz6_new_observ_shap$BaseProbability.sph  <- frz6_new_observ_explain$baseline
frz6_new_observ_shap$FinalProbability.sph <- rowSums(frz6_new_observ_shap)
# Show predicted probabilities as sanity check
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
# Starting point for VUS analysis after writing results previously
##########
# Threshold as determined via Youdens Index on test set, see above
threshold <- 0.286
# Load VUS predictions back in
all_vus_sorted <- read.csv(file=vus_pred_loc)
# Plot 'most benign' and 'most pathogenic' predictions
plotrows <- c(1,2,3,4,5,11217,11218,11219,11220,11221)
for(plotrow in plotrows)
{
  row <- all_vus_sorted[plotrow,]
  p <- shapDecisionPlot(row, threshold)
  pdf_plot_loc <- paste(rootDir, "img", paste0(row$gene, "_", row$delta_aaSeq, ".pdf"), sep="/")
  png_plot_loc <- paste(rootDir, "img", paste0(row$gene, "_", row$delta_aaSeq, ".png"), sep="/")
  # don't save, not using in paper
  #ggsave(filename = pdf_plot_loc, plot = p, device = cairo_pdf, width = 10, height = 4) # height was 6.25
  #ggsave(filename = png_plot_loc, plot = p, device = "png", width = 9.777778, height = 5.5) # 16:9 as PNG for full screen slides
}
# quick 1 off plot for inspection: shapDecisionPlot(all_vus_sorted[11221,], threshold)
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
  p <- shapDecisionPlot(row, threshold)
  pdf_plot_loc <- paste(rootDir, "img", paste0(row$gene, "_", row$delta_aaSeq, ".pdf"), sep="/")
  # don't save all, select specific row later
  #ggsave(filename = pdf_plot_loc, plot = p, device = cairo_pdf, width = 10, height = 4)
}
# as table
vus_changed_sorted <- vus_changed %>% arrange(LP)
vus_changed_sorted$verdict <- ifelse(vus_changed_sorted$LP >= threshold,"P","B")
vus_changed_sorted$dna <- paste0(vus_changed_sorted$dna_variant_chrom,":",vus_changed_sorted$dna_variant_pos," ",vus_changed_sorted$dna_variant_ref,">",vus_changed_sorted$dna_variant_alt)
vus_changed_sorted[,c("gene","TranscriptID","UniProtID","dna","delta_aaSeq","LP", "verdict","new_classification")]
# save select row only
rowWNT7B <- vus_changed_sorted[9,]
p <- shapDecisionPlot(rowWNT7B, threshold)
pdf_plot_loc <- paste(rootDir, "img", paste0("DAVE1_decision_",rowWNT7B$gene, "_", rowWNT7B$delta_aaSeq, ".pdf"), sep="/")
ggsave(filename = pdf_plot_loc, plot = p, device = cairo_pdf, width = 10, height = 4)
# Prep for joint boxplot with ClinVar
vus_changed_sorted_vkgl_min <- vus_changed_sorted[,c("LP", "verdict","new_classification")]
vus_changed_sorted_vkgl_min$new_classification <- ifelse(vus_changed_sorted_vkgl_min$new_classification == "LP","LP/P","LB/B")
vus_changed_sorted_vkgl_min$source <- "VKGL"
# show folding energy change
vus_changed_sorted[,c("gene","UniProtID","delta_aaSeq","LP", "verdict","new_classification", "delta_total.energy")]


#### Now on ClinVar data
# Find variants that were VUS in VKGL April 2024 but have since received a ClinVar classification
# download from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2025/clinvar_20250923.vcf.gz
clinvar_loc <- paste(rootDir, "data", "clinvar_20250923.vcf.gz", sep="/")
clinvarVCF <- read.vcfR(clinvar_loc)
clinvar <- as.data.frame(clinvarVCF@fix)
vus_changed_clinv <- merge(y = clinvar, x = all_vus_sorted, by.y = c("CHROM", "POS", "REF", "ALT"), by.x = c( "dna_variant_chrom", "dna_variant_pos", "dna_variant_ref", "dna_variant_alt"))
# free up memory
clinvarVCF <- NULL
clinvar <- NULL 
# filter clinvar by quality and classifiction
vus_changed_clinv_1star <- subset(vus_changed_clinv, !grepl("CLNREVSTAT=no_assertion_criteria_provided", INFO))
vus_changed_clinv_LP <- subset(vus_changed_clinv_1star, grepl("CLNSIG=(Likely_pathogenic|Pathogenic)", INFO))
vus_changed_clinv_LB <- subset(vus_changed_clinv_1star, grepl("CLNSIG=(Likely_benign|Benign)", INFO))
vus_changed_clinv_LP$new_classification <- "LP/P"
vus_changed_clinv_LB$new_classification <- "LB/B"
vus_changed_clinv_both <- rbind(vus_changed_clinv_LP, vus_changed_clinv_LB)
plot(as.factor(vus_changed_clinv_both$new_classification), vus_changed_clinv_both$FinalProbability.sph)
table(vus_changed_clinv_both$new_classification)
# Prep for joint boxplot with VKGL
vus_changed_sorted_clinvar_min <- vus_changed_clinv_both[,c("LP","new_classification")]
vus_changed_sorted_clinvar_min$verdict <- ifelse(vus_changed_sorted_clinvar_min$LP >= threshold,"P","B")
vus_changed_sorted_clinvar_min$source <- "ClinVar"
# find with affected ligand top pocket
vus_changed_clinv_both_ligand_aff <- vus_changed_clinv_both %>% arrange(delta_ligand_rank1_sas_points)
vus_changed_clinv_both_ligand_aff[c(1,2,3,652,653,654),c("gene","UniProtID","dna","delta_aaSeq","LP","delta_ligand_rank1_sas_points")]
rowLig <- vus_changed_clinv_both_ligand_aff[1,]
p <- shapDecisionPlot(rowLig, threshold)
ligand_plot_loc <- paste(rootDir, "img", paste0("DAVE1_decision_",rowLig$gene, "_", rowLig$delta_aaSeq, ".pdf"), sep="/")
ggsave(filename = ligand_plot_loc, plot = p, device = cairo_pdf, width = 10, height = 4)
# find with affected protein sites
vus_changed_clinv_both_protS_aff <- vus_changed_clinv_both %>% arrange(delta_ProtS_cumu_bin)
vus_changed_clinv_both_protS_aff[c(1,2,3,652,653,654),c("gene","UniProtID","TranscriptID","dna","delta_aaSeq","LP","delta_ProtS_cumu_bin")]
rowProtS <- vus_changed_clinv_both_protS_aff[2,]
p <- shapDecisionPlot(rowProtS, threshold)
rowProtS_plot_loc <- paste(rootDir, "img", paste0("DAVE1_decision_",rowProtS$gene, "_", rowProtS$delta_aaSeq, ".pdf"), sep="/")
ggsave(filename = rowProtS_plot_loc, plot = p, device = cairo_pdf, width = 10, height = 4)

# joint boxplot
shapRed <- "#FF0C57"
shapBlu <- "#1E88E5"
jointData <- rbind(vus_changed_sorted_vkgl_min, vus_changed_sorted_clinvar_min)
jointData$new_classification_by_source <- paste0(jointData$new_classification," in ",jointData$source)
colnames(jointData)
p <- ggplot(jointData, aes(x=new_classification_by_source, y=LP, fill=new_classification)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("LB/B" = shapBlu, "LP/P" = shapRed)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0, face= "italic"),
        plot.subtitle=element_text(size=9),
        plot.title.position = "plot",
        plot.caption.position =  "plot",
        legend.position="none",
        plot.tag.position = c(0.2, 0.025),
        plot.tag=element_text(size=9),
        axis.text.x = element_text(size = 9, colour = "black"),
        axis.text.y = element_text(size = 9, colour = "black"),
        axis.line = element_line(colour = "black")
  ) +
  labs(title = "DAVE1 scores of VKGL April 2024 VUS that have since been classified",
     subtitle = "LB/B = (likely) benign, LP/P = (likely) pathogenic, in VKGL release July 2025 or ClinVar 2025-09-23",
     x = "Expert variant classification per source",
     y = "DAVE1 pathogenicity prediction score")
p
jointbox_plot_loc <- paste(rootDir, "img", paste0("vkgl-clinvar-reclass-jointbox.pdf"), sep="/")
ggsave(filename = jointbox_plot_loc, plot = p, device = cairo_pdf, width = 6, height = 4)



# Apply optimal threshold on VKGL VUS that have since receieved a ClinVar classification
cv <- vus_changed_clinv_both # shorten name
tp <- sum(cv[cv$new_classification=="LP/P",'FinalProbability.sph'] >= threshold)
fp <- sum(cv[cv$new_classification=="LB/B",'FinalProbability.sph'] >= threshold)
tn <- sum(cv[cv$new_classification=="LB/B",'FinalProbability.sph'] < threshold)
fn <- sum(cv[cv$new_classification=="LP/P",'FinalProbability.sph'] < threshold)
ppv <- 100 *tp/(tp+fp)
npv <- 100 *tn/(tn+fn)
sens <- tp / (tp + fn)*100
spec <- tn / (tn + fp)*100
cat(paste("in ClinVar data we applied threshold",threshold,", when applied we find", tp, "TP,", fp, "FP,", tn, "TN and", fn, "FN\n"))
cat(paste("this means ", ppv, "PPV,", npv, "NPV,", sens, "sens and", spec, "spec\n"))

# Apply this threshold on VKGL
#cutpointDF <- subset(vus_changed, new_classification == "LB" | new_classification == "LP")
#opt_cut <- cutpointr(cutpointDF, FinalProbability.sph, new_classification, direction = ">=", pos_class = "LP", neg_class = "LB", method = maximize_metric, metric = youden)
#youdenIndex <- opt_cut$optimal_cutpoint # here, 0.278, but we're using ClinVar's youden
tp <- sum(vus_changed[vus_changed$new_classification=="LP",'FinalProbability.sph'] >= threshold)
fp <- sum(vus_changed[vus_changed$new_classification=="LB",'FinalProbability.sph'] >= threshold)
tn <- sum(vus_changed[vus_changed$new_classification=="LB",'FinalProbability.sph'] < threshold)
fn <- sum(vus_changed[vus_changed$new_classification=="LP",'FinalProbability.sph'] < threshold)
ppv <- 100 *tp/(tp+fp)
npv <- 100 *tn/(tn+fn)
sens <- tp / (tp + fn)*100
spec <- tn / (tn + fp)*100
cat(paste("when applied to VKGL we find", tp, "TP,", fp, "FP,", tn, "TN and", fn, "FN\n"))
cat(paste("this means ", ppv, "PPV,", npv, "NPV,", sens, "sens and", spec, "spec\n"))
# --> 4 TP, 2 FP, 5 TN and 1 FN

# Also, in full VUS set, how many would be over and under the threshold?
nrow(all_vus_sorted[all_vus_sorted$FinalProbability.sph >= threshold, ])
nrow(all_vus_sorted[all_vus_sorted$FinalProbability.sph < threshold, ])
