library(scales)
#library(reshape)
library(reshape2)
library(dplyr)
library(ggplot2)
library(Rmisc) # for CI
library(data.table) # for casting multiple variables
library(cutpointr) # optimal PPV calc


rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability"
#rootDir <- "C:/Users/tk_20/git/vkgl-secretome-protein-stability"
imgDir <- paste(rootDir, "img", sep="/")
setwd(imgDir)

# Load the data and assign meaningful row names
freeze1 <- paste(rootDir, "data", "freeze1.csv", sep="/")
results <- read.csv(freeze1)
rownames(results) <- paste0(results$gene, "/", results$UniProtID, ":", results$delta_aaSeq)


#############################################################################
# Colorblind palette, see https://derekogle.com/NCGraphing/resources/colors #
#############################################################################
LBe <- "#009E73"; VUS <- "#999999"; LPa <- "#D55E00"
sec <- "#E69F00"; int <- "#56B4E9"; mem <- "#CC79A7"
chp <- "#F0E442"; unc <- "#0072B2"; blk <- "#000000"


####################
# Data preparation #
####################
# Remove 'CF' i.e. 'Conflicting' interpretations
results <- subset(results,  ann_classificationVKGL != "CF")
# Remove all columns with only 0 values (for WT, mutant and delta: electrostatic.kon, Entropy.Complex, mloop_entropy, partial.covalent.bonds, sloop_entropy, water.bridge)
results <- results[, colSums(results != 0) > 0]
# Select all columns with relevant factors or numerical variables for analysis
# We drop WT information here since the ann_classificationVKGL applies to the mutant/delta
results <- results %>% select(contains(c("ann_classificationVKGL", "ann_proteinIschaperoned", "ann_proteinLocalization", "mutant_", "delta_")))
results <- results %>% select(-contains(c("_aaSeq", "ann_mutant_energy_SD")))
# Scale all numeric variables to a mean of 0 and SD of 1
results_scaled <- results %>% mutate_if(is.numeric, scale)
# Melt variables except for the factors
results_scaled_melt <- reshape2::melt(results_scaled, na.rm = FALSE, id = c("ann_classificationVKGL", "ann_proteinIschaperoned", "ann_proteinLocalization"))


################################################
# Scatterplot to show localization differences #
################################################
seg_tip_len_scale = 0.01

for(select in c("delta_", "mutant_")){
  # Select only one type of variables to clean up plot
  results_scaled_melt_sub <- results_scaled_melt[grep(select, results_scaled_melt$variable),]
  # Remove variant type tag to clean up plot
  results_scaled_melt_sub$variable = gsub(select, "", results_scaled_melt_sub$variable)
  # Aggregate values on classification, localization and variable
  results_var_means_CI_agg_clsf_loc <- aggregate(results_scaled_melt_sub$value, by=list(Classification=results_scaled_melt_sub$ann_classificationVKGL, Localization=results_scaled_melt_sub$ann_proteinLocalization, Variable=results_scaled_melt_sub$variable), FUN=CI)
  # Cast back so that we have Classification as columns for plotting on X/Y axis
  #results_var_means_CI_agg_clsf_loc_cast <- reshape2::dcast(as.data.table(results_var_means_CI_agg_clsf_loc), Localization+Variable~Classification, value.var=c("x.upper", "x.mean", "x.lower"))
  results_var_means_CI_agg_clsf_loc_cast <- data.table::dcast(as.data.table(results_var_means_CI_agg_clsf_loc), Localization+Variable~Classification, value.var=c("x.upper", "x.mean", "x.lower"))
  # Plot and save
  seg_tip_len_x <- (max(results_var_means_CI_agg_clsf_loc_cast$x.mean_LB)-min(results_var_means_CI_agg_clsf_loc_cast$x.mean_LB))*seg_tip_len_scale
  seg_tip_len_y <- (max(results_var_means_CI_agg_clsf_loc_cast$x.mean_LP)-min(results_var_means_CI_agg_clsf_loc_cast$x.mean_LP))*seg_tip_len_scale
  p <- ggplot(results_var_means_CI_agg_clsf_loc_cast, aes(x =  x.mean_LB, y = x.mean_LP, color = Localization, label=Variable)) +
    theme_classic() +
    geom_point() +
    #geom_abline(intercept = 0, slope = 1) +
    geom_segment(aes(x = x.lower_LB, y = x.mean_LP, xend = x.upper_LB, yend = x.mean_LP)) +
    geom_segment(aes(x = x.mean_LB, y = x.lower_LP, xend = x.mean_LB, yend = x.upper_LP)) +
    geom_segment(aes(x = x.mean_LB-seg_tip_len_x, y = x.upper_LP, xend = x.mean_LB+seg_tip_len_x, yend = x.upper_LP)) + # top tip
    geom_segment(aes(x = x.mean_LB-seg_tip_len_x, y = x.lower_LP, xend = x.mean_LB+seg_tip_len_x, yend = x.lower_LP)) + # bottom tip
    geom_segment(aes(x = x.lower_LB, y = x.mean_LP-seg_tip_len_y, xend = x.lower_LB, yend = x.mean_LP+seg_tip_len_y)) + # left tip
    geom_segment(aes(x = x.upper_LB, y = x.mean_LP-seg_tip_len_y, xend = x.upper_LB, yend = x.mean_LP+seg_tip_len_y)) + # right tip
    geom_text(size=2, hjust = "right", vjust="top", nudge_x = -seg_tip_len_x, nudge_y = -seg_tip_len_y, check_overlap = F) +
    scale_color_manual(name="Protein\nlocalization", labels=c("intracellular" = "Intracellular", "membrane" = "Membrane", "secreted" = "Secreted"), values=c("intracellular"=int, "membrane"=mem, "secreted"=sec))
  ggsave(filename=paste0("scatter-loc-",select,".png"), plot=p, width = 8, height = 4.5)
}

# Same but for chaperoned vs unchaperoned
for(select in c("delta_", "mutant_")){
  results_scaled_melt_sub <- results_scaled_melt[grep(select, results_scaled_melt$variable),]
  results_scaled_melt_sub$variable = gsub(select, "", results_scaled_melt_sub$variable)
  results_var_means_CI_agg_clsf_chp <- aggregate(results_scaled_melt_sub$value, by=list(Classification=results_scaled_melt_sub$ann_classificationVKGL, Chaperoned=results_scaled_melt_sub$ann_proteinIschaperoned, Variable=results_scaled_melt_sub$variable), FUN=CI)
  results_var_means_CI_agg_clsf_chp_cast <- data.table::dcast(as.data.table(results_var_means_CI_agg_clsf_chp), Chaperoned+Variable~Classification, value.var=c("x.upper", "x.mean", "x.lower"))
  seg_tip_len_x <- (max(results_var_means_CI_agg_clsf_chp_cast$x.mean_LB)-min(results_var_means_CI_agg_clsf_chp_cast$x.mean_LB))*seg_tip_len_scale
  seg_tip_len_y <- (max(results_var_means_CI_agg_clsf_chp_cast$x.mean_LP)-min(results_var_means_CI_agg_clsf_chp_cast$x.mean_LP))*seg_tip_len_scale
  p <- ggplot(results_var_means_CI_agg_clsf_chp_cast, aes(x =  x.mean_LB, y = x.mean_LP, color = Chaperoned, label=Variable)) +
    theme_classic() +
    geom_point() +
    #geom_abline(intercept = 0, slope = 1) +
    geom_segment(aes(x = x.lower_LB, y = x.mean_LP, xend = x.upper_LB, yend = x.mean_LP)) +
    geom_segment(aes(x = x.mean_LB, y = x.lower_LP, xend = x.mean_LB, yend = x.upper_LP)) +
    geom_segment(aes(x = x.mean_LB-seg_tip_len_x, y = x.upper_LP, xend = x.mean_LB+seg_tip_len_x, yend = x.upper_LP)) + # top tip
    geom_segment(aes(x = x.mean_LB-seg_tip_len_x, y = x.lower_LP, xend = x.mean_LB+seg_tip_len_x, yend = x.lower_LP)) + # bottom tip
    geom_segment(aes(x = x.lower_LB, y = x.mean_LP-seg_tip_len_y, xend = x.lower_LB, yend = x.mean_LP+seg_tip_len_y)) + # left tip
    geom_segment(aes(x = x.upper_LB, y = x.mean_LP-seg_tip_len_y, xend = x.upper_LB, yend = x.mean_LP+seg_tip_len_y)) + # right tip
    geom_text(size=2, hjust = "right", vjust="top", nudge_x = -seg_tip_len_x, nudge_y = -seg_tip_len_y, check_overlap = F) +
    scale_color_manual(name="Protein is\nchaperoned", labels = c("TRUE" = "Yes", "FALSE" = "No"), values = c("TRUE" = chp, "FALSE" = unc))
  ggsave(filename=paste0("scatter-chp-",select,".png"), plot=p, width = 8, height = 4.5)
}

#################################################
# From this, test specific groups?              #
# Means are meaningless since we split          #
# into LP/LB which often 'cancel eachother out' #
#################################################
chap_all <- subset(results, ann_proteinIschaperoned == TRUE)
unch_all <- subset(results, ann_proteinIschaperoned == FALSE)

secr_all <- subset(results, ann_proteinLocalization == "secreted")
memb_all <- subset(results, ann_proteinLocalization == "membrane")
intr_all <- subset(results, ann_proteinLocalization == "intracellular")

wilcox.test(chap_all$delta_total.energy, unch_all$delta_total.energy)
wilcox.test(secr_all$delta_total.energy, memb_all$delta_total.energy)
wilcox.test(secr_all$delta_total.energy, intr_all$delta_total.energy)
wilcox.test(memb_all$delta_total.energy, intr_all$delta_total.energy)

#################################################################################
# Followup on molecular weight delta for LP variants in different localizations #
#################################################################################
results_LP <- subset(results, ann_classificationVKGL == "LP")
ggplot(results_LP, aes(delta_molWeight, colour = ann_proteinLocalization, fill = ann_proteinLocalization)) +
  theme_classic() +
  geom_density(alpha = 0.25) + #, adjust = 0.1
  scale_fill_manual(name="Protein\nlocalization", labels=c("intracellular" = "Intracellular", "membrane" = "Membrane", "secreted" = "Secreted"), values=c("intracellular"=int, "membrane"=mem, "secreted"=sec)) +
  scale_color_manual(name="Protein\nlocalization", labels=c("intracellular" = "Intracellular", "membrane" = "Membrane", "secreted" = "Secreted"), values=c("intracellular"=int, "membrane"=mem, "secreted"=sec))


######################################################################
# Scatterplot total energy split in localization AND chaperonization #
######################################################################
# Select only total energy
results_scaled_melt_sub <- subset(results_scaled_melt, variable == "delta_total.energy" ) #& ann_classificationVKGL != "VUS"
results_scaled_melt_sub_agg <- aggregate(results_scaled_melt_sub$value, by=list(Classification=results_scaled_melt_sub$ann_classificationVKGL, Localization=results_scaled_melt_sub$ann_proteinLocalization, Chaperoned=results_scaled_melt_sub$ann_proteinIschaperoned, Variable=results_scaled_melt_sub$variable), FUN=CI)
results_scaled_melt_sub_cast <- data.table::dcast(as.data.table(results_scaled_melt_sub_agg), Localization+Chaperoned+Variable~Classification, value.var=c("x.upper", "x.mean", "x.lower"))
seg_tip_len_x <- (max(results_scaled_melt_sub_cast$x.mean_LB)-min(results_scaled_melt_sub_cast$x.mean_LB))*seg_tip_len_scale
seg_tip_len_y <- (max(results_scaled_melt_sub_cast$x.mean_LP)-min(results_scaled_melt_sub_cast$x.mean_LP))*seg_tip_len_scale
p <- ggplot(results_scaled_melt_sub_cast, aes(x = x.mean_LB, y = x.mean_LP, color = Localization, shape = Chaperoned, label=Variable)) +
  theme_classic() +
  geom_point(size=5) +
  geom_segment(aes(x = x.lower_LB, y = x.mean_LP, xend = x.upper_LB, yend = x.mean_LP)) +
  geom_segment(aes(x = x.mean_LB, y = x.lower_LP, xend = x.mean_LB, yend = x.upper_LP)) +
  geom_segment(aes(x = x.mean_LB-seg_tip_len_x, y = x.upper_LP, xend = x.mean_LB+seg_tip_len_x, yend = x.upper_LP)) + # top tip
  geom_segment(aes(x = x.mean_LB-seg_tip_len_x, y = x.lower_LP, xend = x.mean_LB+seg_tip_len_x, yend = x.lower_LP)) + # bottom tip
  geom_segment(aes(x = x.lower_LB, y = x.mean_LP-seg_tip_len_y, xend = x.lower_LB, yend = x.mean_LP+seg_tip_len_y)) + # left tip
  geom_segment(aes(x = x.upper_LB, y = x.mean_LP-seg_tip_len_y, xend = x.upper_LB, yend = x.mean_LP+seg_tip_len_y)) + # right tip
  geom_text(size=2, hjust = "right", vjust="top", nudge_x = -seg_tip_len_x, nudge_y = -seg_tip_len_y, check_overlap = F) +
  scale_color_manual(name="Protein\nlocalization", labels=c("intracellular" = "Intracellular", "membrane" = "Membrane", "secreted" = "Secreted"), values=c("intracellular"=int, "membrane"=mem, "secreted"=sec))
ggsave(filename=paste0("scatter-totalenergy-",select,".png"), plot=p, width = 8, height = 4.5)


###############################################
# Group-specific PPV/NPV thresholds and apply #
###############################################
ppvFromCutpoint <- function(data, cutpoint){
  opt_cut <- cutpointr(data, delta_total.energy, ann_classificationVKGL, direction = ">=", pos_class = "LP", neg_class = "LB", method = cutpointr::oc_manual, cutpoint=cutpoint)
  tp <- sum(data[data$ann_classificationVKGL=="LP","delta_total.energy"] >= cutpoint)
  fp <- sum(data[data$ann_classificationVKGL=="LB","delta_total.energy"] >= cutpoint)
  tn <- sum(data[data$ann_classificationVKGL=="LB","delta_total.energy"] < cutpoint)
  fn <- sum(data[data$ann_classificationVKGL=="LP","delta_total.energy"] < cutpoint)
  ppv <- 100 *tp/(tp+fp)
  npv <- 100 *tn/(tn+fn)
  sens <- opt_cut$sensitivity*100
  spec <- opt_cut$specificity*100
  cat(paste("At cutpoint is ",cutpoint," we get PPV ",round(ppv),"%, NPV ",round(npv),"%, sensitivity ",round(sens),"% and specificity ",round(spec),"%.\n",sep=""))
  return(ppv)
}

findCutpointForPPV <- function(data, minRequestedPPV){
  min <- min(data$delta_total.energy)
  max <- max(data$delta_total.energy)
  for(i in seq(from = min, to = max, length.out = 100)){
    foundPPV <- ppvFromCutpoint(data, i)
    if(foundPPV >= minRequestedPPV)
    return(i)
  }
}

secr_chp <- subset(chap_all, ann_proteinLocalization == "secreted" & ann_classificationVKGL != "VUS")
memb_chp <- subset(chap_all, ann_proteinLocalization == "membrane" & ann_classificationVKGL != "VUS")
intr_chp <- subset(chap_all, ann_proteinLocalization == "intracellular" & ann_classificationVKGL != "VUS")

findCutpointForPPV(secr_chp, 90)
findCutpointForPPV(memb_chp, 90)
findCutpointForPPV(intr_chp, 90)




###########################
# Jitter to show all data #
###########################
p <- ggplot(results_scaled_melt, 
  aes(x = value, 
      y = ann_classificationVKGL,
      color = ann_classificationVKGL)) +
  geom_jitter(size = 0.25) +
  geom_vline(xintercept = 0) +
  facet_grid(rows=vars(variable)) + # cols=vars(ann_proteinLocalization) / ann_proteinIschaperoned
  theme_classic() +
  theme(strip.text.y = element_text(angle = 0), legend.position = "none") +
  scale_colour_manual(name = "Classification", values = c("LB" = "green","LP" = "red", "VUS"="gray"))
ggsave(filename="all-variables-facet.png", plot=p, width = 8, height = 20) #normally 8x4.5
# Check disulfide outliers
highDD <- results_scaled[results_scaled$delta_disulfide > 7,]
table(highDD$ann_classificationVKGL)
paste(rownames(highDD), highDD$ann_classificationVKGL)

