library(ggplot2)
library(plotly)
library(reshape2)
library(grid)
#library(tempR)
library(colorspace)
library(ggrepel)

# Set up locations and working dir
rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability" # root directory that contains README.md, data/, img/, out/, src/, etc.
imgDir <- paste(rootDir, "img", sep="/")
setwd(imgDir)

# Load data
freeze5 <- read.csv(paste(rootDir, "data", "freeze5_predictions.csv", sep="/"))
rownames(freeze5) <- paste0(freeze5$gene, "/", freeze5$UniProtID, ":", freeze5$delta_aaSeq)

# Feature labels and descriptions
feat <- read.csv(paste(rootDir, "data", "features.csv", sep="/"))
feat$Feature <- paste0(feat$Feature, ".sph")
  
# Data slicing
fr5sub <- freeze5[!is.na(freeze5$delta_DNAs_cumu_bin),] # GeoNet not done yet, select where we have data
fr5sub <- subset(fr5sub, MLSplit == "Test") # Train, Test, Unknown - combine with ann_classificationVKGL == LB, LP, VUS, CF

######################
# SHAP Decision Plot #
######################

# Function to to format labels
formatDelta <- function(input) {
  if (is.na(input)) { return("")}
  else if (input > 0) { return(paste("increased by", input))}
  else if (input < 0) { return(paste("decreased by", input))}
  else { return("is unchanged")}
}

# Function to make the SHAP decision plot
# selectRow: data frame row to plot
# featureContribThreshold: terms below will be summed together as one term
shapDecisionPlot <- function(row, featureContribThreshold){
  #row <- freeze5[which(rownames(freeze5) == "GJB1/P08034:CA173Y"),] # DEBUG
  #featureContribThreshold <- 0.01 # DEBUG
  rowSPH <- row[, grepl(".sph$", names(row))] # all rows with a feature SHAP probability heuristic
  rowSPH$FinalProbability # Sanity check pt.1: this cumulative P value should match pt.2 later
  rowSPH <- rowSPH[, !grepl("FinalProbability", names(rowSPH))] # Remove cumulative P from data, we won't use it further
  rowSPHmelt <- reshape2::melt(rowSPH, na.rm = FALSE, id.vars = integer())
  nAboveFeatContribThr <- sum(abs(rowSPHmelt$value) >= featureContribThreshold)
  other <- sum(rowSPHmelt[rev(order(abs(rowSPHmelt$value)))[(nAboveFeatContribThr+1):dim(rowSPHmelt)[1]],]$value) # Order by absolute value and sum the bottom contributors
  rowSPHmelt <- if(nAboveFeatContribThr==0) { rowSPHmelt[0, ] } else { rowSPHmelt[rev(order(abs(rowSPHmelt$value)))[1:nAboveFeatContribThr],] }# Select only the top contributors
  if (!is.na(other)){rowSPHmelt <- rbind(rowSPHmelt, data.frame(variable="OTHER.sph", value=other))} # Add 'other' unless it was NA (i.e. at a threshold of 0)
  rowSPHmelt <- rowSPHmelt[order(abs(rowSPHmelt$value)),] # Re-order by value now that the set is complete
  row.names(rowSPHmelt) <- NULL # Clear out row names (i.e. row indices), resetting them to 1..n
  rowSPHmelt$idx <- as.numeric(row.names(rowSPHmelt)) # Save the row indices so we can always restore this order
  sum(rowSPHmelt$value) # Sanity check pt.2: sum of P values here should match above (ignoring rounding errors)
  rowSPHmelt$prevValue <- c(rowSPHmelt$value[-1], NA) # Add 'previous value' to help plotting (oriented bottom-up)
  
  # Merges with other data
  rowSPHmelt <- merge(rowSPHmelt, feat, by.x = "variable", by.y = "Feature", all.x = T) # Merge with label data
  originalVars <- sub("\\.sph$", "", grep("BaseProbability|OTHER", rowSPHmelt$variable, value = TRUE, invert = TRUE)) # Reconstruct original variable names of current variables except Base and OTHER
  deltaValues <- reshape2::melt(row[originalVars], na.rm = FALSE, id.vars = integer()) # Get original deltas for these variables
  if(ncol(deltaValues) > 0){
    deltaValues$value <- round(deltaValues$value, digits = 6) # Round to prevent values like 1.3000000000109
    deltaValues$variable <- paste0(deltaValues$variable, ".sph") # Equalize variable names for for merging back
    rowSPHmelt <- merge(rowSPHmelt, deltaValues, by.x = "variable", by.y = "variable", all.x = TRUE) # Merge back
    names(rowSPHmelt)[names(rowSPHmelt) == 'value.y'] <- 'delta' # Rename for clarity
    names(rowSPHmelt)[names(rowSPHmelt) == 'value.x'] <- 'value' # Rename back to original
  }
  
  # Restore order and make labels
  rowSPHmelt <- rowSPHmelt[order(rowSPHmelt$idx), ] # Re-order by index since merging swaps things around
  rowSPHmelt$NamePlusEffect <- paste(rowSPHmelt$Name, sapply(rowSPHmelt$delta, formatDelta)) # Create enhanced labels
  rowSPHmelt$idx <- factor(rowSPHmelt$idx, levels = rowSPHmelt$idx, labels = rowSPHmelt$NamePlusEffect) # Assign levels and labels to the indices
  
  # Make raster grob with SHAP colors
  lightenFactor <- 0.33
  nrGradientRows <- nAboveFeatContribThr+1 # Vertical, 1 lane per variable
  nrGradientCols <- 100 # Horizontal, nice to be smooth
  shapRed <- "#FF0C57"
  shapBlu <- "#1E88E5"
  shapRedSoft <- lighten(shapRed, amount = lightenFactor)
  shapBluSoft <- lighten(shapBlu, amount = lightenFactor)
  shapPalette <-colorRampPalette(colors=c(shapBluSoft,shapRedSoft))(nrGradientCols)
  colorMat <- matrix(NA_character_, nrow = nrGradientRows, ncol = nrGradientCols)
  for(i in seq_len(nrGradientCols)) {
    ramp <- colorRampPalette(c("white", shapPalette[i]))(nrGradientRows)
    colorMat[, i] <- ramp
  }
  g <- rasterGrob(colorMat, width = unit(1, "npc"), height = unit(1, "npc"), interpolate = FALSE) # y=0.1, x=0.5,
  
  ## Plot
  rowSPHmelt$valuecs <- rev(cumsum(rev(rowSPHmelt$value)))
  row$ann_classificationVKGL[row$ann_classificationVKGL == "LB"] <- "LB/B"
  row$ann_classificationVKGL[row$ann_classificationVKGL == "LP"] <- "LP/P"
  if(nrow(rowSPHmelt) < 40) { x_text_size <- 10 } else{ x_text_size <- 5 } # NOTE: This text resizing 'cutoff' may change if the number of features in the model changes!
  p <- ggplot(rowSPHmelt, aes(x = idx, y = valuecs, color = value > 0)) +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=0, ymax=1) +
    geom_line(aes(group = 1, color=(value > 0)), linewidth=2, alpha=0.5) +
    geom_vline(xintercept=seq(0, nrow(rowSPHmelt), by=1), linetype="dashed", linewidth=0.2) +
    geom_point() +
    scale_color_manual(labels = c("TRUE" = "More pathogenic", "FALSE" = "More benign"), values = c("TRUE" = shapRed, "FALSE" = shapBlu), name = "Impact") +
    labs(title = paste("SHAP decision plot for ", row$delta_aaSeq, " in gene ", row$gene, " (",row$UniProtID,"), ", row$dna_variant_assembly, " ",
                       row$dna_variant_chrom, ":", row$dna_variant_pos, row$dna_variant_ref, ">", row$dna_variant_alt, ", VKGL April 2024: ", row$ann_classificationVKGL, sep=""),
         subtitle = "SHAP values do not correlate with feature values, but instead capture feature contributions based on the interactions among all features uniquely for this prediction",
         x = paste0("Features affecting probability by >",(featureContribThreshold*100),"%\nin descending order of impact"),
         y = paste("Cumulative probability heuristic for SHAP values\nFinal probability for being pathogenic is", round(sum(rowSPHmelt$value), digits=3))) +
    theme_bw() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.caption = element_text(hjust = 0, face= "italic"),
          plot.subtitle=element_text(size=9),
          plot.title.position = "plot",
          plot.caption.position =  "plot",
          legend.position="none",
          axis.text.x = element_text(size = 9, colour = "black"),
          axis.text.y = element_text(size = x_text_size, colour = "black"), # because of coord flip!
          axis.line = element_line(colour = "black")
          ) + coord_flip()
  ggsave(paste(paste0(row$gene, "_", row$delta_aaSeq, "_thr",featureContribThreshold,".pdf")), plot = p, device = "pdf", width = 10, height = 6.25)
}

# Example function calls
#
# Get top/bottom effects from data slice, retrieve from full set:
# order(fr5sub$FinalProbability.sph, decreasing = T)[1:10]
# fr5sub[x,] -> row name

LP_ExIdx <- which(rownames(freeze5) == "GJB1/P08034:CA201Y")   # LP from Test set
LB_ExIdx <- which(rownames(freeze5) == "CSF2RB/P32927:VA870I") # LB from Test set
VUS_lowP_ExIdx <- which(rownames(freeze5) == "ADAMTSL4/Q6UY14:LA78F") # VUS, low P
VUS_highP_ExIdx <- which(rownames(freeze5) == "GJB1/P08034:CA173Y") # VUS, high P

shapDecisionPlot(freeze5[LP_ExIdx,], 0.00)
shapDecisionPlot(freeze5[LB_ExIdx,], 0.00)
shapDecisionPlot(freeze5[VUS_lowP_ExIdx,], 0.0)
shapDecisionPlot(freeze5[VUS_highP_ExIdx,], 0.0)

shapDecisionPlot(freeze5[LP_ExIdx,], 0.01)
shapDecisionPlot(freeze5[LB_ExIdx,], 0.01)
shapDecisionPlot(freeze5[VUS_lowP_ExIdx,], 0.01)
shapDecisionPlot(freeze5[VUS_highP_ExIdx,], 0.01)



###############
# Other stuff #
###############

# Density
LBe <- "forestgreen"; VUS <- "blue"; LPa <- "red"
p <- ggplot(fr5sub, aes(FinalProbability.sph, colour = ann_classificationVKGL, fill = ann_classificationVKGL)) +
  theme_classic() +
  geom_density(alpha = 0.25,  adjust = 0.3) +
  scale_fill_manual(name="Clsf", labels=c("LB" = "LB", "LP" = "LP", "VUS" = "VUS"), values=c("LB"=LBe, "LP"=LPa, "VUS"=VUS)) +
  scale_color_manual(name="Clsf", labels=c("LB" = "LB", "LP" = "LP", "VUS" = "VUS"), values=c("LB"=LBe, "LP"=LPa, "VUS"=VUS))
ggplotly(p)
ggsave(filename=paste0("train.png"), plot=p, width = 8, height = 4.5)


# VUS reclassif
probCutoff <- 0.5 
predSlice <- subset(fr5sub,  MLSplit == "Test" & FinalProbability.sph > probCutoff )
t <- table(predSlice$ann_classificationVKGL)
t
ppv <- t["LP"]/ (t["LB"] + t["LP"]) 
ppv
subset(predSlice, ann_classificationVKGL=="LB")
vus_to_LP <- subset(fr5sub, ann_classificationVKGL == "VUS" & FinalProbability.sph > probCutoff )
dim(vus_to_LP)
rownames(vus_to_LP)


# compare with AM
p <- ggplot(fr5sub, aes(x = FinalProbability.sph, y = ann_am_pathogenicity, color = ann_classificationVKGL)) +
  theme_classic() +
  geom_point(alpha = 0.25) +
  scale_color_manual(name="Clsf", labels=c("LB" = "LB", "LP" = "LP", "VUS" = "VUS"), values=c("LB"=LBe, "LP"=LPa, "VUS"=VUS))
ggplotly(p)

