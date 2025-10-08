library(colorspace)
library(grid)

# Feature labels and descriptions
feat <- read.csv(paste(rootDir, "data", "12features.csv", sep="/"))
feat$Feature <- paste0(feat$Feature, ".sph")

# Function to to format labels
formatLabels <- function(labelDF) {
  labels <- list()
  for(i in 1:nrow(labelDF)) {
    labelRow <- labelDF[i,]
    if (is.na(labelRow$delta)) {    labels <- c(labels, paste(labelRow["Name"])) }
    else if (labelRow$delta == 0){ labels <- c(labels, paste(labelRow["Name"], "is unaffected", sep=" ")) }
    else if (labelRow$delta > 0) { labels <- c(labels, paste(labelRow["Name"], "increased by", labelRow$delta, unitPlural(labelRow["Unit"], labelRow$delta), sep=" ")) }
    else if (labelRow$delta < 0) { labels <- c(labels, paste(labelRow["Name"], "decreased by", labelRow$delta, unitPlural(labelRow["Unit"], labelRow$delta), sep=" ")) }
    else {                labels <- c(labels, "ERROR: UNDEFINED STATE") }
  }
  return(labels)
}

unitPlural <- function(unit, delta) {
  if(delta==1){return(paste0(unit))}
  if(delta%%1==0){return(paste0(unit,"s"))}
  else{ return(unit)}
}

verdict <- function(finalProb, thrs)
{
  if(finalProb >= thrs){ return("pathogenic") }
  else{ return("benign") }
}

# Function to make the SHAP decision plot
# selectRow: data frame row to plot
shapDecisionPlot <- function(row, thrs){
  featureContribThreshold <- 0 #hardcoded to show all seperately, but can be increased to group contributions
  rowSPH <- row[, grepl(".sph$", names(row))] # all rows with a feature SHAP probability heuristic
  rowSPH$FinalProbability # Sanity check pt.1: this cumulative P value should match pt.2 later
  rowSPH <- rowSPH[, !grepl("FinalProbability", names(rowSPH))] # Remove cumulative P from data, we won't use it further
  rowSPHmelt <- reshape2::melt(rowSPH, na.rm = FALSE, id.vars = integer()) # melt data for further steps
  nAboveFeatContribThr <- sum(abs(rowSPHmelt$value) >= featureContribThreshold) # ascertain nr of features above contribution threshold
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
  rowSPHmelt$NamePlusEffect <- formatLabels(rowSPHmelt) # Create enhanced labels
  rowSPHmelt$idx <- factor(rowSPHmelt$idx, levels = rowSPHmelt$idx, labels = rowSPHmelt$NamePlusEffect) # Assign levels and labels to the indices
  
  # Make raster grob with SHAP colors
  lightenFactor <- 0.33
  nrGradientRows <- nAboveFeatContribThr # Vertical, 1 lane per variable
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
  #row$ann_classificationVKGL[row$ann_classificationVKGL == "LB"] <- "LB/B"
  #row$ann_classificationVKGL[row$ann_classificationVKGL == "LP"] <- "LP/P"
  if(nrow(rowSPHmelt) < 40) { x_text_size <- 10 } else{ x_text_size <- 5 } # NOTE: This text resizing 'cutoff' may change if the number of features in the model changes!
  p <- ggplot(rowSPHmelt, aes(x = idx, y = valuecs, color = value > 0)) +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=0, ymax=1) +
    geom_line(aes(group = 1, color=(value > 0)), linewidth=1, alpha=1, arrow = arrow(angle = 30, length = unit(0.2,"cm"), ends = "first", type = "open")) +
    geom_vline(xintercept=seq(0, nrow(rowSPHmelt), by=1), linetype="dashed", linewidth=0.2) +
    #geom_point() +
    scale_color_manual(labels = c("TRUE" = "More pathogenic", "FALSE" = "More benign"), values = c("TRUE" = shapRed, "FALSE" = shapBlu), name = "Impact") +
    labs(title = paste("DAVE1 decision plot for ", row$delta_aaSeq, " in gene ", row$gene, " (",row$UniProtID,", ",row$TranscriptID,"), ", row$dna_variant_assembly, " ",
                       row$dna_variant_chrom, ":", row$dna_variant_pos, row$dna_variant_ref, ">", row$dna_variant_alt, sep=""),
         subtitle = "SHAP values do not correlate with feature values, but instead capture feature contributions based on the interactions among all features uniquely for this prediction",
         tag = "* = if any, one or more   ** = if any",
         x = paste0("Feature contribution to probability of\nbeing pathogenic, in descending order"),
         y = paste0("Cumulative probability heuristic for SHAP values\nFinal probability for being pathogenic is ", round(sum(rowSPHmelt$value), digits=3), ", variant is estimated to be ", verdict(sum(rowSPHmelt$value), thrs))) +
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
          axis.text.y = element_text(size = x_text_size, colour = "black"), # because of coord flip!
          axis.line = element_line(colour = "black")
    ) + coord_flip()
  return(p)
}
