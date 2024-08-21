library(scales)
#library(reshape)
library(reshape2)
library(dplyr)
library(ggplot2)
library(Rmisc) # for CI
library(data.table) # for casting multiple variables


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
# Select only 'mutant_' variables to clean up plot
results_scaled_melt_sub <- results_scaled_melt[grep("mutant_", results_scaled_melt$variable),]
# Remove 'mutant_' tag to clean up plot
results_scaled_melt_sub$variable = gsub("mutant_", "", results_scaled_melt_sub$variable)
# Aggregate values on classification, localization and variable
results_var_means_CI_agg_clsf_loc <- aggregate(results_scaled_melt_sub$value, by=list(Classification=results_scaled_melt_sub$ann_classificationVKGL, Localization=results_scaled_melt_sub$ann_proteinLocalization, Variable=results_scaled_melt_sub$variable), FUN=CI)
# Cast back so that we have Classification as columns for plotting on X/Y axis
#results_var_means_CI_agg_clsf_loc_cast <- reshape2::dcast(as.data.table(results_var_means_CI_agg_clsf_loc), Localization+Variable~Classification, value.var=c("x.upper", "x.mean", "x.lower"))
results_var_means_CI_agg_clsf_loc_cast <- data.table::dcast(as.data.table(results_var_means_CI_agg_clsf_loc), Localization+Variable~Classification, value.var=c("x.upper", "x.mean", "x.lower"))
head(results_var_means_CI_agg_clsf_loc_cast)
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
  geom_text(size=2, hjust = "left", vjust="top", nudge_x = 0.01, nudge_y = -0.01, check_overlap = TRUE) +
  scale_color_manual(name="Protein localization", labels=c("intracellular" = "Intracellular", "membrane" = "Membrane", "secreted" = "Secreted"), values=c("intracellular"=int, "membrane"=mem, "secreted"=sec)) +
  xlim(-1.1, 1.1) +
  ylim(-1.1, 1.1)
ggsave(filename="loc-scatter-mut.png", plot=p, width = 8, height = 4.5)

# Same but for deltas, has no limits
results_scaled_melt_sub <- results_scaled_melt[grep("delta_", results_scaled_melt$variable),]
results_scaled_melt_sub$variable = gsub("delta_", "", results_scaled_melt_sub$variable)
results_var_means_CI_agg_clsf_loc <- aggregate(results_scaled_melt_sub$value, by=list(Classification=results_scaled_melt_sub$ann_classificationVKGL, Localization=results_scaled_melt_sub$ann_proteinLocalization, Variable=results_scaled_melt_sub$variable), FUN=CI)
results_var_means_CI_agg_clsf_loc_cast <- data.table::dcast(as.data.table(results_var_means_CI_agg_clsf_loc), Localization+Variable~Classification, value.var=c("x.upper", "x.mean", "x.lower"))
head(results_var_means_CI_agg_clsf_loc_cast)
seg_tip_len_x <- (max(results_var_means_CI_agg_clsf_loc_cast$x.mean_LB)-min(results_var_means_CI_agg_clsf_loc_cast$x.mean_LB))*seg_tip_len_scale
seg_tip_len_y <- (max(results_var_means_CI_agg_clsf_loc_cast$x.mean_LP)-min(results_var_means_CI_agg_clsf_loc_cast$x.mean_LP))*seg_tip_len_scale
p <- ggplot(results_var_means_CI_agg_clsf_loc_cast, aes(x =  x.mean_LB, y = x.mean_LP, color = Localization, label=Variable)) +
  theme_classic() +
  geom_point() +
  geom_segment(aes(x = x.lower_LB, y = x.mean_LP, xend = x.upper_LB, yend = x.mean_LP)) +
  geom_segment(aes(x = x.mean_LB, y = x.lower_LP, xend = x.mean_LB, yend = x.upper_LP)) +
  geom_segment(aes(x = x.mean_LB-seg_tip_len_x, y = x.upper_LP, xend = x.mean_LB+seg_tip_len_x, yend = x.upper_LP)) + # top tip
  geom_segment(aes(x = x.mean_LB-seg_tip_len_x, y = x.lower_LP, xend = x.mean_LB+seg_tip_len_x, yend = x.lower_LP)) + # bottom tip
  geom_segment(aes(x = x.lower_LB, y = x.mean_LP-seg_tip_len_y, xend = x.lower_LB, yend = x.mean_LP+seg_tip_len_y)) + # left tip
  geom_segment(aes(x = x.upper_LB, y = x.mean_LP-seg_tip_len_y, xend = x.upper_LB, yend = x.mean_LP+seg_tip_len_y)) + # right tip
  geom_text(size=2, hjust = "left", vjust="top", nudge_x = 0.01, nudge_y = -0.01, check_overlap = TRUE) +
  scale_color_manual(name="Protein localization", labels=c("intracellular" = "Intracellular", "membrane" = "Membrane", "secreted" = "Secreted"), values=c("intracellular"=int, "membrane"=mem, "secreted"=sec))
ggsave(filename="loc-scatter-delta.png", plot=p, width = 8, height = 4.5)

# Same but for chaperoned vs unchaperoned, select mutants/deltas
for(select in c("delta_", "mutant_")){
  results_scaled_melt_sub <- results_scaled_melt[grep(select, results_scaled_melt$variable),]
  results_scaled_melt_sub$variable = gsub(select, "", results_scaled_melt_sub$variable)
  results_var_means_CI_agg_clsf_chp <- aggregate(results_scaled_melt_sub$value, by=list(Classification=results_scaled_melt_sub$ann_classificationVKGL, Chaperoned=results_scaled_melt_sub$ann_proteinIschaperoned, Variable=results_scaled_melt_sub$variable), FUN=CI)
  results_var_means_CI_agg_clsf_chp_cast <- data.table::dcast(as.data.table(results_var_means_CI_agg_clsf_chp), Chaperoned+Variable~Classification, value.var=c("x.upper", "x.mean", "x.lower"))
  head(results_var_means_CI_agg_clsf_chp_cast)
  seg_tip_len_x <- (max(results_var_means_CI_agg_clsf_loc_cast$x.mean_LB)-min(results_var_means_CI_agg_clsf_loc_cast$x.mean_LB))*seg_tip_len_scale
  seg_tip_len_y <- (max(results_var_means_CI_agg_clsf_loc_cast$x.mean_LP)-min(results_var_means_CI_agg_clsf_loc_cast$x.mean_LP))*seg_tip_len_scale
  p <- ggplot(results_var_means_CI_agg_clsf_chp_cast, aes(x =  x.mean_LB, y = x.mean_LP, color = Chaperoned, label=Variable)) +
    theme_classic() +
    geom_point() +
    #geom_abline(intercept = 0, slope = 1) +
    geom_segment(aes(x = x.lower_LB, y = x.mean_LP, xend = x.upper_LB, yend = x.mean_LP)) +
    geom_segment(aes(x = x.mean_LB, y = x.lower_LP, xend = x.mean_LB, yend = x.upper_LP)) +
    geom_segment(aes(x = x.mean_LB-seg_tip_len_x, y = x.upper_LP, xend = x.mean_LB+seg_tip_len_x, yend = x.upper_LP)) + # top tip
    geom_segment(aes(x = x.mean_LB-seg_tip_len_x, y = x.lower_LP, xend = x.mean_LB+seg_tip_len_x, yend = x.lower_LP)) + # bottom tip
    geom_segment(aes(x = x.lower_LB, y = x.mean_LP-seg_tip_len, xend = x.lower_LB, yend = x.mean_LP+seg_tip_len_y)) + # left tip
    geom_segment(aes(x = x.upper_LB, y = x.mean_LP-seg_tip_len_y, xend = x.upper_LB, yend = x.mean_LP+seg_tip_len_y)) + # right tip
    geom_text(size=2, hjust = "left", vjust="top", nudge_x = 0.01, nudge_y = -0.01, check_overlap = TRUE) +
    scale_color_manual(values = c("TRUE" = chp, "FALSE" = unc))
  ggsave(filename=paste0("chp-scatter-",select,".png"), plot=p, width = 8, height = 4.5)
}


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

