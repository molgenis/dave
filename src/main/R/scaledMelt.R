library(scales)
#library(reshape)
library(reshape2)
library(dplyr)
library(ggplot2)

rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability"
imgDir <- paste(rootDir, "img", sep="/")
setwd(imgDir)

freeze1 <- paste(rootDir, "data", "freeze1.csv", sep="/")
results <- read.csv(freeze1)
results <- subset(results,  ann_classificationVKGL != "CF")
rownames(results) <- paste0(results$gene, "/", results$UniProtID, ":", results$delta_aaSeq)
# Select all relevant columns for plot
results <- results %>% select(contains(c("ann_classificationVKGL", "ann_proteinIschaperoned", "ann_proteinLocalization", "delta_", "mutant_")))
results <- results %>% select(-contains(c("_aaSeq", "WT_", "ann_mutant_energy_SD")))
# Remove null columns
results <- results[, colSums(results != 0) > 0]
# Scale all values in each column to a mean of 0 and SD of 1
results_scaled <- results %>% mutate_if(is.numeric, scale)
# Melt variables except for the 'factors'
results_scaled_melt <- reshape2::melt(results_scaled, na.rm = FALSE, id = c("ann_classificationVKGL", "ann_proteinIschaperoned", "ann_proteinLocalization"))
# Jitter plot
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


highDD <- results_scaled[results_scaled$delta_disulfide > 7,]
table(highDD$ann_classificationVKGL)
paste(rownames(highDD), highDD$ann_classificationVKGL)
