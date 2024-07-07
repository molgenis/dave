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


######################################################################
# Wild-type folding energy vs. mass with localization and chaperoned #
######################################################################
ggplot(resGeneColl, aes(x=wtDG, y=mwDa, color=protType, shape=chaperoned, label=gene)) +
  theme_classic() +
  geom_point() +
  geom_text(size = 2, hjust=-0.1, vjust=-0.3, check_overlap = TRUE)+
  scale_shape_manual(name="Chaperoned", values=c(1,3), labels=c("FALSE"="No", "TRUE"="Yes")) +
  scale_color_manual(name="Protein localization", labels=c("intracellular" = "Intracellular", "membrane" = "Membrane", "secreted" = "Secreted"), values=c("intracellular"=int, "membrane"=mem, "secreted"=sec)) +
  scale_y_continuous(labels = label_comma()) +
  scale_x_continuous(labels = label_comma()) +
  xlab("Gibbs free energy change of wild-type protein folding (in kcal/mol)") +
  ylab("Protein molecular weight (in Daltons)")
ggsave("protein_wtDG_vs_mwDa_scatterplot.png", width = 8, height = 4.5)


###########################################
# Groupwise comparisons with significance #
###########################################
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


####################################################
# Panel of things
####################################################
axisTextSize <- 7
axisTitleSize <- 9
titleSize <- 9
memb_intr_unc <- rbind(memb_unc, intr_unc)
p1 <- ggplot(data = memb_intr_unc, aes(x = protType, y = wtDG, fill = protType)) +
  theme_classic() +
  geom_violin() +
  ggtitle(paste0("Adj. p-value = ", format.pval(stats[1,]$adjPvalue, 3), " ", stats[1,]$adjPstars)) +
  scale_fill_manual(values = c("intracellular" = int, "membrane" = mem)) +
  scale_x_discrete(name="Protein localization, unchaperoned", labels=c(
  "membrane" = paste0("Membrane",
         "\nmedian = ", prettyNum(median(memb_intr_unc[which(memb_intr_unc$protType=="membrane"),]$mwDa), big.mark = ","),
         "\nmean = ", prettyNum(round(mean(memb_intr_unc[which(memb_intr_unc$protType=="membrane"),]$mwDa)), big.mark = ",")
  ),
  "intracellular" = paste0("Intracellular",
         "\nmedian = ", prettyNum(median(memb_intr_unc[which(memb_intr_unc$protType=="intracellular"),]$mwDa), big.mark = ","),
         "\nmean = ", prettyNum(round(mean(memb_intr_unc[which(memb_intr_unc$protType=="intracellular"),]$mwDa)), big.mark = ",")
  )))+
  theme(axis.text=element_text(size=axisTextSize, color="black"), axis.title=element_text(size=axisTitleSize), plot.title = element_text(size = titleSize), legend.position = "none") +
  scale_y_continuous(labels = label_comma()) +
  ylab("Wild-type ΔG")

secr_memb <- rbind(secr_all, memb_all)
p2 <- ggplot(data = secr_memb, aes(x = protType, y = mwDa, fill = protType)) +
  theme_classic() +
  geom_violin() +
  ggtitle(paste0("Adj. p-value = ", format.pval(stats[2,]$adjPvalue, 3), " ", stats[2,]$adjPstars)) +
  scale_fill_manual(values = c("secreted" = sec, "membrane" = mem)) +
  scale_x_discrete(name="Protein localization", labels=c(
    "membrane" = paste0("Membrane",
           "\nmedian = ", prettyNum(median(secr_memb[which(secr_memb$protType=="membrane"),]$mwDa), big.mark = ","),
           "\nmean = ", prettyNum(round(mean(secr_memb[which(secr_memb$protType=="membrane"),]$mwDa)), big.mark = ",")
    ),
    "secreted" = paste0("Secreted",
           "\nmedian = ", prettyNum(median(secr_memb[which(secr_memb$protType=="secreted"),]$mwDa), big.mark = ","),
           "\nmean = ", prettyNum(round(mean(secr_memb[which(secr_memb$protType=="secreted"),]$mwDa)), big.mark = ",")
    )))+
  theme(axis.text=element_text(size=axisTextSize, color="black"), axis.title=element_text(size=axisTitleSize), plot.title = element_text(size = titleSize), legend.position = "none") +
  scale_y_continuous(labels = label_comma()) +
  ylab("Mass")

memb_intr <- rbind(memb_all, intr_all)
p3 <- ggplot(data = memb_intr, aes(x = protType, y = wtDG, fill = protType)) +
  theme_classic() +
  geom_violin() +
  ggtitle(paste0("Adj. p-value = ", format.pval(stats[3,]$adjPvalue, 3), " ", stats[3,]$adjPstars)) +
  scale_fill_manual(values = c("intracellular" = int, "membrane" = mem)) +
  scale_x_discrete(name="Protein localization", labels=c(
    "membrane" = paste0("Membrane",
           "\nmedian = ", prettyNum(median(memb_intr[which(memb_intr$protType=="membrane"),]$wtDG), big.mark = ","),
           "\nmean = ", prettyNum(round(mean(memb_intr[which(memb_intr$protType=="membrane"),]$wtDG)), big.mark = ",")
    ),
    "intracellular" = paste0("Intracellular",
           "\nmedian = ", prettyNum(median(memb_intr[which(memb_intr$protType=="intracellular"),]$wtDG), big.mark = ","),
           "\nmean = ", prettyNum(round(mean(memb_intr[which(memb_intr$protType=="intracellular"),]$wtDG)), big.mark = ",")
    )))+
  theme(axis.text=element_text(size=axisTextSize, color="black"), axis.title=element_text(size=axisTitleSize), plot.title = element_text(size = titleSize), legend.position = "none") +
  scale_y_continuous(labels = label_comma()) +
  ylab("Wild-type ΔG")

secr_memb_chp <- rbind(secr_chp, memb_chp)
p4 <- ggplot(data = secr_memb_chp, aes(x = protType, y = mwDa, fill = protType)) +
  theme_classic() +
  geom_violin() +
  ggtitle(paste0("Adj. p-value = ", format.pval(stats[4,]$adjPvalue, 3), " ", stats[4,]$adjPstars)) +
  scale_fill_manual(values = c("secreted" = sec, "membrane" = mem)) +
  scale_x_discrete(name="Protein localization, chaperoned", labels=c(
    "membrane" = paste0("Membrane",
                      "\nmedian = ", prettyNum(median(secr_memb_chp[which(secr_memb_chp$protType=="membrane"),]$mwDa), big.mark = ","),
                      "\nmean = ", prettyNum(round(mean(secr_memb_chp[which(secr_memb_chp$protType=="membrane"),]$mwDa)), big.mark = ",")
    ),
    "secreted" = paste0("Secreted",
                           "\nmedian = ", prettyNum(median(secr_memb_chp[which(secr_memb_chp$protType=="secreted"),]$mwDa), big.mark = ","),
                           "\nmean = ", prettyNum(round(mean(secr_memb_chp[which(secr_memb_chp$protType=="secreted"),]$mwDa)), big.mark = ",")
    )))+
  theme(axis.text=element_text(size=axisTextSize, color="black"), axis.title=element_text(size=axisTitleSize), plot.title = element_text(size = titleSize), legend.position = "none") +
  scale_y_continuous(labels = label_comma()) +
  ylab("Mass")

secr_intr <- rbind(secr_all, intr_all)
p5 <- ggplot(data = secr_intr, aes(x = protType, y = mwDa, fill = protType)) +
  theme_classic() +
  geom_violin() +
  ggtitle(paste0("Adj. p-value = ", format.pval(stats[5,]$adjPvalue, 3), " ", stats[5,]$adjPstars)) +
  scale_fill_manual(values = c("secreted" = sec, "intracellular" = mem)) +
  scale_x_discrete(name="Protein localization", labels=c(
    "intracellular" = paste0("Intracellular",
                        "\nmedian = ", prettyNum(median(secr_intr[which(secr_intr$protType=="intracellular"),]$mwDa), big.mark = ","),
                        "\nmean = ", prettyNum(round(mean(secr_intr[which(secr_intr$protType=="intracellular"),]$mwDa)), big.mark = ",")
    ),
    "secreted" = paste0("Secreted",
                        "\nmedian = ", prettyNum(median(secr_intr[which(secr_intr$protType=="secreted"),]$mwDa), big.mark = ","),
                        "\nmean = ", prettyNum(round(mean(secr_intr[which(secr_intr$protType=="secreted"),]$mwDa)), big.mark = ",")
    )))+
  theme(axis.text=element_text(size=axisTextSize, color="black"), axis.title=element_text(size=axisTitleSize), plot.title = element_text(size = titleSize), legend.position = "none") +
  scale_y_continuous(labels = label_comma()) +
  ylab("Mass")

secr_intr_unc <- rbind(secr_unc, intr_unc)
p6 <- ggplot(data = secr_intr_unc, aes(x = protType, y = wtDG, fill = protType)) +
  theme_classic() +
  geom_violin() +
  ggtitle(paste0("Adj. p-value = ", format.pval(stats[6,]$adjPvalue, 3), " ", stats[6,]$adjPstars)) +
  scale_fill_manual(values = c("secreted" = sec, "intracellular" = int)) +
  scale_x_discrete(name="Protein localization, unchaperoned", labels=c(
    "intracellular" = paste0("Intracellular",
                             "\nmedian = ", prettyNum(median(secr_intr_unc[which(secr_intr_unc$protType=="intracellular"),]$wtDG), big.mark = ","),
                             "\nmean = ", prettyNum(round(mean(secr_intr_unc[which(secr_intr_unc$protType=="intracellular"),]$wtDG)), big.mark = ",")
    ),
    "secreted" = paste0("Secreted",
                        "\nmedian = ", prettyNum(median(secr_intr_unc[which(secr_intr_unc$protType=="secreted"),]$wtDG), big.mark = ","),
                        "\nmean = ", prettyNum(round(mean(secr_intr_unc[which(secr_intr_unc$protType=="secreted"),]$wtDG)), big.mark = ",")
    )))+
  theme(axis.text=element_text(size=axisTextSize, color="black"), axis.title=element_text(size=axisTitleSize), plot.title = element_text(size = titleSize), legend.position = "none") +
  scale_y_continuous(labels = label_comma()) +
  ylab("Wild-type ΔG")

patchwork <- (p1 + p2 + p3) / (p4 + p5 + p6)
patchwork + plot_annotation(tag_levels = 'A', 
                            title = 'Top 6 most significant differences in protein properties',
                            subtitle = 'Axis labels: \'Mass\' = protein molecular mass in Daltons. \'Wild-type ΔG\' = Gibbs free energy change of wild-type protein folding in kcal/mol.',
                            theme = theme(plot.subtitle = element_text(size = 8)))
ggsave("patchwork.png", width = 8, height = 4.5)


