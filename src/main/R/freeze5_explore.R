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
fr5sub <- subset(freeze5, MLSplit != "Train") # Train, Test, Unknown / LB, LP, VUS

######################
# SHAP Decision Plot #
######################

# Data prep
topImportantFeatures <- 65 # 66 max. Alternative: implement automatic '% contribution' cutoff
row <- fr5sub[1,] # 134 / 217. see top effects with order(fr5sub$FinalProbability.sph, decreasing = T)[1:10]
rowPS <- row[, grepl(".sph$", names(row))]
rowPS$FinalProbability # sanity check: this should match ...
rowPS <- rowPS[, !grepl("FinalProbability", names(rowPS))]
rowPSmelt <- reshape2::melt(rowPS, na.rm = FALSE)
other <- sum(rowPSmelt[rev(order(abs(rowPSmelt$value)))[(topImportantFeatures+1):dim(rowPSmelt)[1]],]$value)
rowPSmelt <- rowPSmelt[rev(order(abs(rowPSmelt$value)))[1:topImportantFeatures],]
rowPSmelt <- rbind(rowPSmelt, data.frame(variable="OTHER.sph", value=other))
rowPSmelt <- rowPSmelt[order(abs(rowPSmelt$value)),]
row.names(rowPSmelt) <- NULL
rowPSmelt$idx <- as.numeric(row.names(rowPSmelt))
rowPSmelt$prevValue <- c(rowPSmelt$value[-1], NA)
sum(rowPSmelt$value) # ... this! (ignoring rounding errors)
rowPSmelt <- merge(rowPSmelt, feat, by.x = "variable", by.y = "Feature", all.x = T) # merge with label data
rowPSmelt <- rowPSmelt[order(rowPSmelt$idx), ] # re-order after merge
rowPSmelt$idx <- factor(rowPSmelt$idx, levels = rowPSmelt$idx, labels = rowPSmelt$Name)

# Predefined raster grob with SHAP colors
lightenFactor <- 0.33
nrGradientRows <- topImportantFeatures+1
nrGradientCols <- 7
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
rowPSmelt$valuecs <- rev(cumsum(rev(rowPSmelt$value)))
rowPSmelt$valuecs[rowPSmelt$valuecs > 1] <- 1 # rounding errors can lead to >1, fix here
ggplot(rowPSmelt, aes(x = idx, y = valuecs, color = value > 0)) +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) + 
  geom_line(aes(group = 1, color=(value > 0)), linewidth=2, alpha=0.5) +
  geom_point() +
  scale_color_manual(labels = c("TRUE" = "More pathogenic", "FALSE" = "More benign"), values = c("TRUE" = shapRed, "FALSE" = shapBlu), name = "Impact") +
  labs(title = "SHAP Decision Plot for a Single Observation", x = "x", y = "y") +
  theme_minimal() + theme(legend.position="none") + ylim(0,1) + coord_flip() 




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

