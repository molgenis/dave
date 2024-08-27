library(ggplot2)   # for plotting
library(ggpubr)    # for boxplots with pvalues
library(ggbeeswarm)# for adding dots to plots
library(dplyr)     # to remove duplicate rows
library(scales)    # for big values with commas in plots
library(ggrepel)   # alternative for geom_text: geom_text_repel
library(patchwork) # multiple plots in one
library(gtools)    # pvals to stars


##################################
# Directories and data locations #
##################################
rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability"
freeze1 <- paste(rootDir, "data", "freeze1.csv", sep="/")
results <- read.csv(freeze1)
imgDir <- paste(rootDir, "img", sep="/")
setwd(imgDir)


#############################################################################
# Colorblind palette, see https://derekogle.com/NCGraphing/resources/colors #
#############################################################################
LBe <- "#009E73"; VUS <- "#999999"; LPa <- "#D55E00"
sec <- "#E69F00"; int <- "#56B4E9"; mem <- "#CC79A7"
chp <- "#F0E442"; unc <- "#0072B2"; blk <- "#000000"


#########################################################
# Groupwise comparisons LP/LB, chaperoned, localization #
#########################################################

# Primary groups - subsetting
LB_all <- subset(results, classificationVKGL == "LB")
LP_all <- subset(results, classificationVKGL == "LP")

chap_all <- subset(results, chaperoned == TRUE)
unch_all <- subset(results, chaperoned == FALSE)

secr_all <- subset(results, protType == "secreted")
memb_all <- subset(results, protType == "membrane")
intr_all <- subset(results, protType == "intracellular")

# Secondary groups - subsetting
chap_LB <- subset(LB_all, chaperoned == TRUE)
unch_LB <- subset(LB_all, chaperoned == FALSE)
chap_LP <- subset(LP_all, chaperoned == TRUE)
unch_LP <- subset(LP_all, chaperoned == FALSE)

chap_secr <- subset(secr_all, chaperoned == TRUE)
unch_secr <- subset(secr_all, chaperoned == FALSE)
chap_memb <- subset(memb_all, chaperoned == TRUE)
unch_memb <- subset(memb_all, chaperoned == FALSE)
chap_intr <- subset(intr_all, chaperoned == TRUE)
unch_intr <- subset(intr_all, chaperoned == FALSE)

secr_LB <- subset(LB_all, protType == "secreted")
secr_LP <- subset(LP_all, protType == "secreted")
memb_LB <- subset(LB_all, protType == "membrane")
memb_LP <- subset(LP_all, protType == "membrane")
intr_LB <- subset(LB_all, protType == "intracellular")
intr_LP <- subset(LP_all, protType == "intracellular")

# Tertiary groups - subsetting
secr_chp_LB <- subset(secr_LB, chaperoned == TRUE)
secr_unc_LB <- subset(secr_LB, chaperoned == FALSE)
secr_chp_LP <- subset(secr_LP, chaperoned == TRUE)
secr_unc_LP <- subset(secr_LP, chaperoned == FALSE)

memb_chp_LB <- subset(memb_LB, chaperoned == TRUE)
memb_unc_LB <- subset(memb_LB, chaperoned == FALSE)
memb_chp_LP <- subset(memb_LP, chaperoned == TRUE)
memb_unc_LP <- subset(memb_LP, chaperoned == FALSE)

intr_chp_LB <- subset(intr_LB, chaperoned == TRUE)
intr_unc_LB <- subset(intr_LB, chaperoned == FALSE)
intr_chp_LP <- subset(intr_LP, chaperoned == TRUE)
intr_unc_LP <- subset(intr_LP, chaperoned == FALSE)

# Primary groups - testing
LB_VS_LP <- wilcox.test(LB_all$total.energy, LP_all$total.energy)
chap_VS_unch <- wilcox.test(chap_all$total.energy, unch_all$total.energy)
secr_VS_memb <- wilcox.test(secr_all$total.energy, memb_all$total.energy)
secr_VS_intr <- wilcox.test(secr_all$total.energy, intr_all$total.energy)
memb_VS_intr <- wilcox.test(memb_all$total.energy, intr_all$total.energy)

# Secondary groups - testing
chap_LB_vs_unch_LB <- wilcox.test(chap_LB$total.energy, unch_LB$total.energy)
chap_LP_vs_unch_LP <- wilcox.test(chap_LP$total.energy, unch_LP$total.energy)
chap_LB_vs_chap_LP <- wilcox.test(chap_LB$total.energy, chap_LP$total.energy)
unch_LB_vs_unch_LP <- wilcox.test(unch_LB$total.energy, unch_LP$total.energy)

chap_secr_VS_unch_secr <- wilcox.test(chap_secr$total.energy, unch_secr$total.energy)
chap_secr_VS_chap_memb <- wilcox.test(chap_secr$total.energy, chap_memb$total.energy)
chap_secr_VS_chap_intr <- wilcox.test(chap_secr$total.energy, chap_intr$total.energy)
unch_secr_VS_unch_memb <- wilcox.test(unch_secr$total.energy, unch_memb$total.energy)
unch_secr_VS_unch_intr <- wilcox.test(unch_secr$total.energy, unch_intr$total.energy)
chap_memb_VS_unch_memb <- wilcox.test(chap_memb$total.energy, unch_memb$total.energy)
chap_memb_VS_chap_intr <- wilcox.test(chap_memb$total.energy, chap_intr$total.energy)
unch_memb_VS_unch_intr <- wilcox.test(unch_memb$total.energy, unch_intr$total.energy)
chap_intr_VS_unch_intr <- wilcox.test(chap_intr$total.energy, unch_intr$total.energy)

secr_LB_VS_secr_LP <- wilcox.test(secr_LB$total.energy, secr_LP$total.energy)
secr_LB_VS_memb_LB <- wilcox.test(secr_LB$total.energy, memb_LB$total.energy)
secr_LB_VS_intr_LB <- wilcox.test(secr_LB$total.energy, intr_LB$total.energy)
secr_LP_VS_memb_LP <- wilcox.test(secr_LP$total.energy, memb_LP$total.energy)
secr_LP_VS_intr_LP <- wilcox.test(secr_LP$total.energy, intr_LP$total.energy)
memb_LB_VS_memb_LP <- wilcox.test(memb_LB$total.energy, memb_LP$total.energy)
memb_LB_VS_intr_LB <- wilcox.test(memb_LB$total.energy, intr_LB$total.energy)
memb_LP_VS_intr_LP <- wilcox.test(memb_LP$total.energy, intr_LP$total.energy)
intr_LB_VS_intr_LP <- wilcox.test(intr_LB$total.energy, intr_LP$total.energy)

# Tertiary groups - testing
secr_chp_LB_VS_secr_chp_LP <- wilcox.test(secr_chp_LB$total.energy, secr_chp_LP$total.energy)
secr_chp_LB_VS_secr_unc_LB <- wilcox.test(secr_chp_LB$total.energy, secr_unc_LB$total.energy)
secr_chp_LB_VS_memb_chp_LB <- wilcox.test(secr_chp_LB$total.energy, memb_chp_LB$total.energy)
secr_chp_LB_VS_intr_chp_LB <- wilcox.test(secr_chp_LB$total.energy, intr_chp_LB$total.energy)
secr_unc_LB_VS_secr_unc_LP <- wilcox.test(secr_unc_LB$total.energy, secr_unc_LP$total.energy)
secr_unc_LB_VS_memb_unc_LB <- wilcox.test(secr_unc_LB$total.energy, memb_unc_LB$total.energy)
secr_unc_LB_VS_intr_unc_LB <- wilcox.test(secr_unc_LB$total.energy, intr_unc_LB$total.energy)
secr_chp_LP_VS_secr_unc_LP <- wilcox.test(secr_chp_LP$total.energy, secr_unc_LP$total.energy)
secr_chp_LP_VS_memb_chp_LP <- wilcox.test(secr_chp_LP$total.energy, memb_chp_LP$total.energy)
secr_chp_LP_VS_intr_chp_LP <- wilcox.test(secr_chp_LP$total.energy, intr_chp_LP$total.energy)
secr_unc_LP_VS_memb_unc_LP <- wilcox.test(secr_unc_LP$total.energy, memb_unc_LP$total.energy)
secr_unc_LP_VS_intr_unc_LP <- wilcox.test(secr_unc_LP$total.energy, intr_unc_LP$total.energy)
memb_chp_LB_VS_memb_chp_LP <- wilcox.test(memb_chp_LB$total.energy, memb_chp_LP$total.energy)
memb_chp_LB_VS_memb_unc_LB <- wilcox.test(memb_chp_LB$total.energy, memb_unc_LB$total.energy)
memb_chp_LB_VS_intr_chp_LB <- wilcox.test(memb_chp_LB$total.energy, intr_chp_LB$total.energy)
memb_unc_LB_VS_memb_unc_LP <- wilcox.test(memb_unc_LB$total.energy, memb_unc_LP$total.energy)
memb_unc_LB_VS_intr_unc_LB <- wilcox.test(memb_unc_LB$total.energy, intr_unc_LB$total.energy)
memb_chp_LP_VS_memb_unc_LP <- wilcox.test(memb_chp_LP$total.energy, memb_unc_LP$total.energy)
memb_chp_LP_VS_intr_chp_LP <- wilcox.test(memb_chp_LP$total.energy, intr_chp_LP$total.energy)
memb_unc_LP_VS_intr_unc_LP <- wilcox.test(memb_unc_LP$total.energy, intr_unc_LP$total.energy)
intr_chp_LB_VS_intr_chp_LP <- wilcox.test(intr_chp_LB$total.energy, intr_chp_LP$total.energy)
intr_chp_LB_VS_intr_unc_LB <- wilcox.test(intr_chp_LB$total.energy, intr_unc_LB$total.energy)
intr_unc_LB_VS_intr_unc_LP <- wilcox.test(intr_unc_LB$total.energy, intr_unc_LP$total.energy)
intr_chp_LP_VS_intr_unc_LP  <- wilcox.test(intr_chp_LP$total.energy, intr_unc_LP $total.energy)

# Primary groups - result bundling
stats <- data.frame(Comparison=character(), Variable=character(), pValue=numeric(), stars=character())
stats <- rbind(stats, data.frame(Comparison="LB vs LP", Variable="ΔΔG", pValue=LB_VS_LP$p.value, pStars=stars.pval(LB_VS_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="Chaperoned vs unchaperoned", Variable="ΔΔG", pValue=chap_VS_unch$p.value, pStars=stars.pval(chap_VS_unch$p.value)))
stats <- rbind(stats, data.frame(Comparison="Secreted vs membrane", Variable="ΔΔG", pValue=secr_VS_memb$p.value, pStars=stars.pval(secr_VS_memb$p.value)))
stats <- rbind(stats, data.frame(Comparison="Secreted vs intracellular", Variable="ΔΔG", pValue=secr_VS_intr$p.value, pStars=stars.pval(secr_VS_intr$p.value)))
stats <- rbind(stats, data.frame(Comparison="Membrane vs intracellular", Variable="ΔΔG", pValue=memb_VS_intr$p.value, pStars=stars.pval(memb_VS_intr$p.value)))

# Secondary groups - result bundling
stats <- rbind(stats, data.frame(Comparison="chap_LB_vs_unch_LB", Variable="ΔΔG", pValue=chap_LB_vs_unch_LB$p.value, pStars=stars.pval(chap_LB_vs_unch_LB$p.value)))
stats <- rbind(stats, data.frame(Comparison="chap_LP_vs_unch_LP", Variable="ΔΔG", pValue=chap_LP_vs_unch_LP$p.value, pStars=stars.pval(chap_LP_vs_unch_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="chap_LB_vs_chap_LP", Variable="ΔΔG", pValue=chap_LB_vs_chap_LP$p.value, pStars=stars.pval(chap_LB_vs_chap_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="unch_LB_vs_unch_LP", Variable="ΔΔG", pValue=unch_LB_vs_unch_LP$p.value, pStars=stars.pval(unch_LB_vs_unch_LP$p.value)))

stats <- rbind(stats, data.frame(Comparison="chap_secr_VS_unch_secr", Variable="ΔΔG", pValue=chap_secr_VS_unch_secr$p.value, pStars=stars.pval(chap_secr_VS_unch_secr$p.value)))
stats <- rbind(stats, data.frame(Comparison="chap_secr_VS_chap_memb", Variable="ΔΔG", pValue=chap_secr_VS_chap_memb$p.value, pStars=stars.pval(chap_secr_VS_chap_memb$p.value)))
stats <- rbind(stats, data.frame(Comparison="chap_secr_VS_chap_intr", Variable="ΔΔG", pValue=chap_secr_VS_chap_intr$p.value, pStars=stars.pval(chap_secr_VS_chap_intr$p.value)))
stats <- rbind(stats, data.frame(Comparison="unch_secr_VS_unch_memb", Variable="ΔΔG", pValue=unch_secr_VS_unch_memb$p.value, pStars=stars.pval(unch_secr_VS_unch_memb$p.value)))
stats <- rbind(stats, data.frame(Comparison="unch_secr_VS_unch_intr", Variable="ΔΔG", pValue=unch_secr_VS_unch_intr$p.value, pStars=stars.pval(unch_secr_VS_unch_intr$p.value)))
stats <- rbind(stats, data.frame(Comparison="chap_memb_VS_unch_memb", Variable="ΔΔG", pValue=chap_memb_VS_unch_memb$p.value, pStars=stars.pval(chap_memb_VS_unch_memb$p.value)))
stats <- rbind(stats, data.frame(Comparison="chap_memb_VS_chap_intr", Variable="ΔΔG", pValue=chap_memb_VS_chap_intr$p.value, pStars=stars.pval(chap_memb_VS_chap_intr$p.value)))
stats <- rbind(stats, data.frame(Comparison="unch_memb_VS_unch_intr", Variable="ΔΔG", pValue=unch_memb_VS_unch_intr$p.value, pStars=stars.pval(unch_memb_VS_unch_intr$p.value)))
stats <- rbind(stats, data.frame(Comparison="chap_intr_VS_unch_intr", Variable="ΔΔG", pValue=chap_intr_VS_unch_intr$p.value, pStars=stars.pval(chap_intr_VS_unch_intr$p.value)))

stats <- rbind(stats, data.frame(Comparison="secr_LB_VS_secr_LP", Variable="ΔΔG", pValue=secr_LB_VS_secr_LP$p.value, pStars=stars.pval(secr_LB_VS_secr_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="secr_LB_VS_memb_LB", Variable="ΔΔG", pValue=secr_LB_VS_memb_LB$p.value, pStars=stars.pval(secr_LB_VS_memb_LB$p.value)))
stats <- rbind(stats, data.frame(Comparison="secr_LB_VS_intr_LB", Variable="ΔΔG", pValue=secr_LB_VS_intr_LB$p.value, pStars=stars.pval(secr_LB_VS_intr_LB$p.value)))
stats <- rbind(stats, data.frame(Comparison="secr_LP_VS_memb_LP", Variable="ΔΔG", pValue=secr_LP_VS_memb_LP$p.value, pStars=stars.pval(secr_LP_VS_memb_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="secr_LP_VS_intr_LP", Variable="ΔΔG", pValue=secr_LP_VS_intr_LP$p.value, pStars=stars.pval(secr_LP_VS_intr_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="memb_LB_VS_memb_LP", Variable="ΔΔG", pValue=memb_LB_VS_memb_LP$p.value, pStars=stars.pval(memb_LB_VS_memb_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="memb_LB_VS_intr_LB", Variable="ΔΔG", pValue=memb_LB_VS_intr_LB$p.value, pStars=stars.pval(memb_LB_VS_intr_LB$p.value)))
stats <- rbind(stats, data.frame(Comparison="memb_LP_VS_intr_LP", Variable="ΔΔG", pValue=memb_LP_VS_intr_LP$p.value, pStars=stars.pval(memb_LP_VS_intr_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="intr_LB_VS_intr_LP", Variable="ΔΔG", pValue=intr_LB_VS_intr_LP$p.value, pStars=stars.pval(intr_LB_VS_intr_LP$p.value)))

# Tertiary groups - result bundling
stats <- rbind(stats, data.frame(Comparison="secr_chp_LB_VS_secr_chp_LP", Variable="ΔΔG", pValue=secr_chp_LB_VS_secr_chp_LP$p.value, pStars=stars.pval(secr_chp_LB_VS_secr_chp_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="secr_chp_LB_VS_secr_unc_LB", Variable="ΔΔG", pValue=secr_chp_LB_VS_secr_unc_LB$p.value, pStars=stars.pval(secr_chp_LB_VS_secr_unc_LB$p.value)))
stats <- rbind(stats, data.frame(Comparison="secr_chp_LB_VS_memb_chp_LB", Variable="ΔΔG", pValue=secr_chp_LB_VS_memb_chp_LB$p.value, pStars=stars.pval(secr_chp_LB_VS_memb_chp_LB$p.value)))
stats <- rbind(stats, data.frame(Comparison="secr_chp_LB_VS_intr_chp_LB", Variable="ΔΔG", pValue=secr_chp_LB_VS_intr_chp_LB$p.value, pStars=stars.pval(secr_chp_LB_VS_intr_chp_LB$p.value)))
stats <- rbind(stats, data.frame(Comparison="secr_unc_LB_VS_secr_unc_LP", Variable="ΔΔG", pValue=secr_unc_LB_VS_secr_unc_LP$p.value, pStars=stars.pval(secr_unc_LB_VS_secr_unc_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="secr_unc_LB_VS_memb_unc_LB", Variable="ΔΔG", pValue=secr_unc_LB_VS_memb_unc_LB$p.value, pStars=stars.pval(secr_unc_LB_VS_memb_unc_LB$p.value)))
stats <- rbind(stats, data.frame(Comparison="secr_unc_LB_VS_intr_unc_LB", Variable="ΔΔG", pValue=secr_unc_LB_VS_intr_unc_LB$p.value, pStars=stars.pval(secr_unc_LB_VS_intr_unc_LB$p.value)))
stats <- rbind(stats, data.frame(Comparison="secr_chp_LP_VS_secr_unc_LP", Variable="ΔΔG", pValue=secr_chp_LP_VS_secr_unc_LP$p.value, pStars=stars.pval(secr_chp_LP_VS_secr_unc_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="secr_chp_LP_VS_memb_chp_LP", Variable="ΔΔG", pValue=secr_chp_LP_VS_memb_chp_LP$p.value, pStars=stars.pval(secr_chp_LP_VS_memb_chp_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="secr_chp_LP_VS_intr_chp_LP", Variable="ΔΔG", pValue=secr_chp_LP_VS_intr_chp_LP$p.value, pStars=stars.pval(secr_chp_LP_VS_intr_chp_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="secr_unc_LP_VS_memb_unc_LP", Variable="ΔΔG", pValue=secr_unc_LP_VS_memb_unc_LP$p.value, pStars=stars.pval(secr_unc_LP_VS_memb_unc_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="secr_unc_LP_VS_intr_unc_LP", Variable="ΔΔG", pValue=secr_unc_LP_VS_intr_unc_LP$p.value, pStars=stars.pval(secr_unc_LP_VS_intr_unc_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="memb_chp_LB_VS_memb_chp_LP", Variable="ΔΔG", pValue=memb_chp_LB_VS_memb_chp_LP$p.value, pStars=stars.pval(memb_chp_LB_VS_memb_chp_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="memb_chp_LB_VS_memb_unc_LB", Variable="ΔΔG", pValue=memb_chp_LB_VS_memb_unc_LB$p.value, pStars=stars.pval(memb_chp_LB_VS_memb_unc_LB$p.value)))
stats <- rbind(stats, data.frame(Comparison="memb_chp_LB_VS_intr_chp_LB", Variable="ΔΔG", pValue=memb_chp_LB_VS_intr_chp_LB$p.value, pStars=stars.pval(memb_chp_LB_VS_intr_chp_LB$p.value)))
stats <- rbind(stats, data.frame(Comparison="memb_unc_LB_VS_memb_unc_LP", Variable="ΔΔG", pValue=memb_unc_LB_VS_memb_unc_LP$p.value, pStars=stars.pval(memb_unc_LB_VS_memb_unc_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="memb_unc_LB_VS_intr_unc_LB", Variable="ΔΔG", pValue=memb_unc_LB_VS_intr_unc_LB$p.value, pStars=stars.pval(memb_unc_LB_VS_intr_unc_LB$p.value)))
stats <- rbind(stats, data.frame(Comparison="memb_chp_LP_VS_memb_unc_LP", Variable="ΔΔG", pValue=memb_chp_LP_VS_memb_unc_LP$p.value, pStars=stars.pval(memb_chp_LP_VS_memb_unc_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="memb_chp_LP_VS_intr_chp_LP", Variable="ΔΔG", pValue=memb_chp_LP_VS_intr_chp_LP$p.value, pStars=stars.pval(memb_chp_LP_VS_intr_chp_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="memb_unc_LP_VS_intr_unc_LP", Variable="ΔΔG", pValue=memb_unc_LP_VS_intr_unc_LP$p.value, pStars=stars.pval(memb_unc_LP_VS_intr_unc_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="intr_chp_LB_VS_intr_chp_LP", Variable="ΔΔG", pValue=intr_chp_LB_VS_intr_chp_LP$p.value, pStars=stars.pval(intr_chp_LB_VS_intr_chp_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="intr_chp_LB_VS_intr_unc_LB", Variable="ΔΔG", pValue=intr_chp_LB_VS_intr_unc_LB$p.value, pStars=stars.pval(intr_chp_LB_VS_intr_unc_LB$p.value)))
stats <- rbind(stats, data.frame(Comparison="intr_unc_LB_VS_intr_unc_LP", Variable="ΔΔG", pValue=intr_unc_LB_VS_intr_unc_LP$p.value, pStars=stars.pval(intr_unc_LB_VS_intr_unc_LP$p.value)))
stats <- rbind(stats, data.frame(Comparison="intr_chp_LP_VS_intr_unc_LP ", Variable="ΔΔG", pValue=intr_chp_LP_VS_intr_unc_LP $p.value, pStars=stars.pval(intr_chp_LP_VS_intr_unc_LP $p.value)))

stats$adjPvalue <- ifelse(stats$pValue*nrow(stats) < 1, stats$pValue*nrow(stats), 1)
stats$adjPstars <- stars.pval(stats$adjPvalue)
stats <- stats[order(stats$pValue, decreasing = FALSE), ]   
stats


####
# Panel
###
axisTextSize <- 7
axisTitleSize <- 9
titleSize <- 9

makeBenignVsPathogenicPlot <- function(theData, title){
  plot <- ggplot(data = theData, aes(x = classificationVKGL, y = total.energy, fill = classificationVKGL)) +
    theme_classic() +
    geom_violin() +
    ggtitle(paste0("Adj. p-value = ", format.pval(stats[1,]$adjPvalue, 3), " ", stats[1,]$adjPstars)) +
    scale_fill_manual(values = c("LP" = LPa, "LB" = LBe)) +
    scale_x_discrete(name=title, labels=c(
      "LP" = paste0("Pathogenic or likely pathogenic",
                    "\nmedian = ", round(median(theData[which(theData$classificationVKGL=="LP"),]$total.energy), 2),
                    "\nmean = ", round(mean(theData[which(theData$classificationVKGL=="LP"),]$total.energy), 2)
      ),
      "LB" = paste0("Benign or likely benign",
                    "\nmedian = ", round(median(theData[which(theData$classificationVKGL=="LB"),]$total.energy), 2),
                    "\nmean = ", round(mean(theData[which(theData$classificationVKGL=="LB"),]$total.energy), 2)
      )))+
    theme(axis.text=element_text(size=axisTextSize, color="black"), axis.title=element_text(size=axisTitleSize), plot.title = element_text(size = titleSize), legend.position = "none") +
    scale_y_continuous(labels = label_comma()) +
    ylab("ΔΔG")
  return(plot)
}

p1 <- makeBenignVsPathogenicPlot(rbind(LB_all, LP_all), "Consensus variant classification VKGL release April 2024")
p2 <- makeBenignVsPathogenicPlot(rbind(chap_LB, chap_LP), "Chaperoned proteins, variant classification")
p3 <- makeBenignVsPathogenicPlot(rbind(unch_LB, unch_LP), "Unchaperoned proteins, variant classification")
p4 <- makeBenignVsPathogenicPlot(rbind(secr_LB, secr_LP), "Secreted proteins, variant classification")
p5 <- makeBenignVsPathogenicPlot(rbind(memb_LB, memb_LP), "Membrane proteins, variant classification")
p6 <- makeBenignVsPathogenicPlot(rbind(intr_LB, intr_LP), "Intracellular proteins, variant classification")
p1 / (p2 + p3) / (p4 + p5 + p6)


#####
# Scatterplots
#####
# ann_classificationVKGL, ann_proteinLocalization, ann_proteinIschaperoned, WT_total.energy, ..
lplb <- subset(results, ann_classificationVKGL == "LP" | ann_classificationVKGL == "LB")
ggplot(lplb %>% arrange(match(ann_classificationVKGL, c("LB", "LP"))), aes(y=delta_total.energy, x=WT_total.energy, color=ann_classificationVKGL)) +
  theme_classic() +
  geom_point(size=1) +
  scale_color_manual(name="VKGL classification", values=c("LB" = "green", "LP" = "red")) +
  #scale_color_manual(name="Protein\nlocalization", labels=c("intracellular" = "Intracellular", "membrane" = "Membrane", "secreted" = "Secreted"), values=c("intracellular"=int, "membrane"=mem, "secreted"=sec)) +
  scale_x_continuous(labels = label_comma()) +
  xlab("Gibbs free energy change of wild-type protein folding (ΔG, in kcal/mol)") +
  ylab(paste("Difference in Gibbs free energy change between wild-type and variant protein (ΔΔG, in kcal/mol) ", sep=" "))
# Mol weight?
ggplot(lplb, aes(y=delta_molWeight, x=WT_total.energy, color=ann_proteinLocalization)) +
  theme_classic() +
  geom_point(size=1)

