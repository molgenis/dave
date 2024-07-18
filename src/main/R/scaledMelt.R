library(scales)
library(reshape)
library(reshape2)
library(dplyr)
library(ggplot2)

rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability"
imgDir <- paste(rootDir, "img", sep="/")
setwd(imgDir)

freeze1 <- paste(rootDir, "data", "freeze1.csv", sep="/")
results <- read.csv(freeze1)
# Select all relevant columns for plot
results <- results %>% select(contains(c("ann_classificationVKGL", "ann_proteinIschaperoned", "ann_proteinLocalization", "delta_", "mutant_")))
results <- results %>% select(-contains(c("_aaSeq", "WT_", "ann_mutant_energy_SD")))
# Remove null columns
results <- results[, colSums(results != 0) > 0]
# Scale all values in each column to a mean of 0 and SD of 1
results <- results %>% mutate_if(is.numeric, scale)
# Melt variables except for the 'factors'
results_melt <- reshape2::melt(results, na.rm = FALSE, id = c("ann_classificationVKGL", "ann_proteinIschaperoned", "ann_proteinLocalization"))

# Subset for plot
lplb <- subset(results_melt, ann_classificationVKGL == "LP" | ann_classificationVKGL == "LB")

ggplot(lplb %>% arrange(match(ann_classificationVKGL, c("LB", "LP"))), aes(y=variable, x=value, color=ann_classificationVKGL)) +
  theme_classic() +
  geom_jitter(size=0.1) +
  scale_color_manual(values = c("LB" = "green", "LP" = "red"))
ggsave("all-variables-scaled.png", width = 8, height = 20) #normally 8x4.5

