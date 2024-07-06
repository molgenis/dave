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


##################################################
# Collapse variant-level results into gene-level #
##################################################
resGeneColl <- data.frame(gene=results$gene, mwDa=results$mwDa, wtDG=results$wtDG, transcript=results$transcript, uniprot=results$uniprot, protType=results$protType, chaperoned=results$chaperoned)
resGeneColl <- resGeneColl %>% distinct()
resGeneColl$energyPerKDa <- resGeneColl$wtDG/(resGeneColl$mwDa/1000)


#############################################################################
# Colorblind palette, see https://derekogle.com/NCGraphing/resources/colors #
#############################################################################
LBe <- "#009E73"; VUS <- "#999999"; LPa <- "#D55E00"
sec <- "#E69F00"; int <- "#56B4E9"; mem <- "#CC79A7"
chp <- "#F0E442"; unc <- "#0072B2"; blk <- "#000000"


####
# Stats
####
secr_all <- subset(resGeneColl, protType == "secreted")
secr_chp <- subset(secr_all, chaperoned == TRUE)
secr_unc <- subset(secr_all, chaperoned == FALSE)
memb_all <- subset(resGeneColl, protType == "membrane")
memb_chp <- subset(memb_all, chaperoned == TRUE)
memb_unc <- subset(memb_all, chaperoned == FALSE)
intr_all <- subset(resGeneColl, protType == "intracellular")
intr_chp <- subset(intr_all, chaperoned == TRUE)
intr_unc <- subset(intr_all, chaperoned == FALSE)

secr_all_mwDa_VS_memb_all_mwDa <- wilcox.test(secr_all$mwDa, memb_all$mwDa)
secr_all_wtDG_VS_memb_all_wtDG <- wilcox.test(secr_all$wtDG, memb_all$wtDG)
secr_all_mwDa_VS_intr_all_mwDa <- wilcox.test(secr_all$mwDa, intr_all$mwDa)
secr_all_wtDG_VS_intr_all_wtDG <- wilcox.test(secr_all$wtDG, intr_all$wtDG)
memb_all_mwDa_VS_intr_all_mwDa <- wilcox.test(memb_all$mwDa, intr_all$mwDa)
memb_all_wtDG_VS_intr_all_wtDG <- wilcox.test(memb_all$wtDG, intr_all$wtDG)

secr_chp_mwDa_VS_memb_chp_mwDa <- wilcox.test(secr_chp$mwDa, memb_chp$mwDa)
secr_unc_mwDa_VS_memb_unc_mwDa <- wilcox.test(secr_unc$mwDa, memb_unc$mwDa)
secr_chp_wtDG_VS_memb_chp_wtDG <- wilcox.test(secr_chp$wtDG, memb_chp$wtDG)
secr_unc_wtDG_VS_memb_unc_wtDG <- wilcox.test(secr_unc$wtDG, memb_unc$wtDG)

secr_chp_mwDa_VS_intr_chp_mwDa <- wilcox.test(secr_chp$mwDa, intr_chp$mwDa)
secr_unc_mwDa_VS_intr_unc_mwDa <- wilcox.test(secr_unc$mwDa, intr_unc$mwDa)
secr_chp_wtDG_VS_intr_chp_wtDG <- wilcox.test(secr_chp$wtDG, intr_chp$wtDG)
secr_unc_wtDG_VS_intr_unc_wtDG <- wilcox.test(secr_unc$wtDG, intr_unc$wtDG)

memb_chp_mwDa_VS_intr_chp_mwDa <- wilcox.test(memb_chp$mwDa, intr_chp$mwDa)
memb_unc_mwDa_VS_intr_unc_mwDa <- wilcox.test(memb_unc$mwDa, intr_unc$mwDa)
memb_chp_wtDG_VS_intr_chp_wtDG <- wilcox.test(memb_chp$wtDG, intr_chp$wtDG)
memb_unc_wtDG_VS_intr_unc_wtDG <- wilcox.test(memb_unc$wtDG, intr_unc$wtDG)

stats <- data.frame(Comparison=character(), Variable=character(), pValue=numeric(), stars=character())
stats <- rbind(stats, data.frame(Comparison="Secreted vs membrane", Variable="Molecular mass", pValue=secr_all_mwDa_VS_memb_all_mwDa$p.value, pStars=stars.pval(secr_all_mwDa_VS_memb_all_mwDa$p.value)))
stats <- rbind(stats, data.frame(Comparison="Secreted vs membrane", Variable="Folding energy", pValue=secr_all_wtDG_VS_memb_all_wtDG$p.value, pStars=stars.pval(secr_all_wtDG_VS_memb_all_wtDG$p.value)))
stats <- rbind(stats, data.frame(Comparison="Secreted vs intracellular", Variable="Molecular mass", pValue=secr_all_mwDa_VS_intr_all_mwDa$p.value, pStars=stars.pval(secr_all_mwDa_VS_intr_all_mwDa$p.value)))
stats <- rbind(stats, data.frame(Comparison="Secreted vs intracellular", Variable="Folding energy", pValue=secr_all_wtDG_VS_intr_all_wtDG$p.value, pStars=stars.pval(secr_all_wtDG_VS_intr_all_wtDG$p.value)))
stats <- rbind(stats, data.frame(Comparison="Membrane vs intracellular", Variable="Molecular mass", pValue=memb_all_mwDa_VS_intr_all_mwDa$p.value, pStars=stars.pval(memb_all_mwDa_VS_intr_all_mwDa$p.value)))
stats <- rbind(stats, data.frame(Comparison="Membrane vs intracellular", Variable="Folding energy", pValue=memb_all_wtDG_VS_intr_all_wtDG$p.value, pStars=stars.pval(memb_all_wtDG_VS_intr_all_wtDG$p.value)))

stats <- rbind(stats, data.frame(Comparison="Secreted chaperoned vs membrane chaperoned", Variable="Molecular mass", pValue=secr_chp_mwDa_VS_memb_chp_mwDa$p.value, pStars=stars.pval(secr_chp_mwDa_VS_memb_chp_mwDa$p.value)))
stats <- rbind(stats, data.frame(Comparison="Secreted unchaperoned vs membrane unchaperoned", Variable="Molecular mass", pValue=secr_unc_mwDa_VS_memb_unc_mwDa$p.value, pStars=stars.pval(secr_unc_mwDa_VS_memb_unc_mwDa$p.value)))
stats <- rbind(stats, data.frame(Comparison="Secreted chaperoned vs membrane chaperoned", Variable="Folding energy", pValue=secr_chp_wtDG_VS_memb_chp_wtDG$p.value, pStars=stars.pval(secr_chp_wtDG_VS_memb_chp_wtDG$p.value)))
stats <- rbind(stats, data.frame(Comparison="Secreted unchaperoned vs membrane unchaperoned", Variable="Folding energy", pValue=secr_unc_wtDG_VS_memb_unc_wtDG$p.value, pStars=stars.pval(secr_unc_wtDG_VS_memb_unc_wtDG$p.value)))

stats <- rbind(stats, data.frame(Comparison="Secreted chaperoned vs intracellular chaperoned", Variable="Molecular mass", pValue=secr_chp_mwDa_VS_intr_chp_mwDa$p.value, pStars=stars.pval(secr_chp_mwDa_VS_intr_chp_mwDa$p.value)))
stats <- rbind(stats, data.frame(Comparison="Secreted unchaperoned vs intracellular unchaperoned", Variable="Molecular mass", pValue=secr_unc_mwDa_VS_intr_unc_mwDa$p.value, pStars=stars.pval(secr_unc_mwDa_VS_intr_unc_mwDa$p.value)))
stats <- rbind(stats, data.frame(Comparison="Secreted chaperoned vs intracellular chaperoned", Variable="Folding energy", pValue=secr_chp_wtDG_VS_intr_chp_wtDG$p.value, pStars=stars.pval(secr_chp_wtDG_VS_intr_chp_wtDG$p.value)))
stats <- rbind(stats, data.frame(Comparison="Secreted unchaperoned vs intracellular unchaperoned", Variable="Folding energy", pValue=secr_unc_wtDG_VS_intr_unc_wtDG$p.value, pStars=stars.pval(secr_unc_wtDG_VS_intr_unc_wtDG$p.value)))

stats <- rbind(stats, data.frame(Comparison="Membrane chaperoned vs intracellular chaperoned", Variable="Molecular mass", pValue=memb_chp_mwDa_VS_intr_chp_mwDa$p.value, pStars=stars.pval(memb_chp_mwDa_VS_intr_chp_mwDa$p.value)))
stats <- rbind(stats, data.frame(Comparison="Membrane unchaperoned vs intracellular unchaperoned", Variable="Molecular mass", pValue=memb_unc_mwDa_VS_intr_unc_mwDa$p.value, pStars=stars.pval(memb_unc_mwDa_VS_intr_unc_mwDa$p.value)))
stats <- rbind(stats, data.frame(Comparison="Membrane chaperoned vs intracellular chaperoned", Variable="Folding energy", pValue=memb_chp_wtDG_VS_intr_chp_wtDG$p.value, pStars=stars.pval(memb_chp_wtDG_VS_intr_chp_wtDG$p.value)))
stats <- rbind(stats, data.frame(Comparison="Membrane unchaperoned vs intracellular unchaperoned", Variable="Folding energy", pValue=memb_unc_wtDG_VS_intr_unc_wtDG$p.value, pStars=stars.pval(memb_unc_wtDG_VS_intr_unc_wtDG$p.value)))

stats$adjPvalue <- ifelse(stats$pValue*nrow(stats) < 1, stats$pValue*nrow(stats), 1)
stats$adjPstars <- stars.pval(stats$adjPvalue)
stats <- stats[order(stats$pValue, decreasing = FALSE), ]   
stats
# Top 6:
stats[1:6, ]


