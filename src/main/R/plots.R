library(ggplot2)   # for plotting
library(ggpubr)    # for boxplots with pvalues
library(ggbeeswarm)# for adding dots to plots
library(dplyr)     # to remove duplicate rows
library(scales)    # for big values with commas in plots

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
ex1 <- "#F0E442"; ex2 <- "#0072B2"; ex3 <- "#000000"


####################################################
# Difference in mass between protein localizations #
####################################################
compare <- list( c("intracellular", "membrane"), c("intracellular", "secreted"), c("secreted", "membrane") )
ggviolin(resGeneColl, x = "protType", y = "mwDa", color = "protType", fill= "protType")+ 
  geom_quasirandom(size=0.5)+
  scale_fill_manual(values = c(sec, mem, int)) +
  stat_kruskal_test(label.y = 470000, label = "Kruskal-Wallis p = {p} {p.signif}") +
  geom_pwc(method = "wilcox_test", label = "Wilcoxon adj. p = {p.adj} {p.adj.signif}") +
  scale_x_discrete(name="Protein localization", labels=c(
    paste0("Secreted",
           "\nmedian = ", prettyNum(median(resGeneColl[which(resGeneColl$protType=="secreted"),]$mwDa), big.mark = ","),
           "\nmean = ", prettyNum(round(mean(resGeneColl[which(resGeneColl$protType=="secreted"),]$mwDa)), big.mark = ",")
    ),
    paste0("Membrane",
           "\nmedian = ", prettyNum(median(resGeneColl[which(resGeneColl$protType=="membrane"),]$mwDa), big.mark = ","),
           "\nmean = ", prettyNum(round(mean(resGeneColl[which(resGeneColl$protType=="membrane"),]$mwDa)), big.mark = ",")
    ),
    paste0("Intracellular",
           "\nmedian = ", prettyNum(median(resGeneColl[which(resGeneColl$protType=="intracellular"),]$mwDa), big.mark = ","),
           "\nmean = ", prettyNum(round(mean(resGeneColl[which(resGeneColl$protType=="intracellular"),]$mwDa)), big.mark = ",")
    )
  ))+
  scale_color_manual(values=c(sec, mem, int)) +
  theme(axis.text=element_text(size=10)) +
  scale_y_continuous(labels = label_comma()) +
  ylab("Protein molecular mass (in Daltons)") +
  theme(legend.position = "none")
ggsave("mwDa_difference_test.png", width = 8, height = 4.5)


########################################################################
# Difference in wild-type folding energy between protein localizations #
########################################################################
ggviolin(resGeneColl, x = "protType", y = "wtDG", color = "protType", fill= "protType")+ 
  geom_quasirandom(size=0.01, width=0.23)+
  scale_fill_manual(values = c(sec, mem, int)) +
  stat_kruskal_test(label.y = 7950, label = "Kruskal-Wallis p = {p} {p.signif}") +
  geom_pwc(method = "wilcox_test", label = "Wilcoxon adj. p = {p.adj} {p.adj.signif}") +
  scale_x_discrete(name="Protein localization", labels=c(
    paste0("Secreted",
           "\nmedian = ", round(median(resGeneColl[which(resGeneColl$protType=="secreted"),]$wtDG),2),
           "\nmean = ", round(mean(resGeneColl[which(resGeneColl$protType=="secreted"),]$wtDG),2)
          ),
    paste0("Membrane",
          "\nmedian = ", round(median(resGeneColl[which(resGeneColl$protType=="membrane"),]$wtDG),2),
          "\nmean = ", round(mean(resGeneColl[which(resGeneColl$protType=="membrane"),]$wtDG),2)
          ),
    paste0("Intracellular",
          "\nmedian = ", round(median(resGeneColl[which(resGeneColl$protType=="intracellular"),]$wtDG),2),
          "\nmean = ", round(mean(resGeneColl[which(resGeneColl$protType=="intracellular"),]$wtDG),2)
          )
    ))+
  scale_color_manual(values=c(sec, mem, int)) +
  theme(axis.text=element_text(size=10)) +
  scale_y_continuous(labels = label_comma()) +
  ylab("Gibbs free energy change of wild-type\nprotein folding (in kcal/mol)") +
  theme(legend.position = "none")
ggsave("wtDG_difference_test.png", width = 8, height = 4.5)








kruskal.test(wtDG ~ protType, data = resGeneColl)
kruskal.test(mwDa ~ protType, data = resGeneColl)
ic <- subset(resGeneColl, protType == "intracellular")
mb <- subset(resGeneColl, protType == "membrane")
sc <- subset(resGeneColl, protType == "secreted")
mean(ic$wtDG)
mean(mb$wtDG)
mean(sc$wtDG)
mean(ic$mwDa)
mean(mb$mwDa)
mean(sc$mwDa)
ggplot(resGeneColl, aes(x=wtDG, y=mwDa, color=protType, shape=chaperoned,label=gene)) +
  theme_classic() +
  geom_point() +
  geom_text(size = 3, hjust=-0.1, vjust=-0.1, check_overlap = TRUE)+
  scale_shape_manual(values=c(1,3)) +
  scale_color_manual(values=c("purple", "blue", "black")) +
  scale_y_continuous(labels = label_comma()) +
  xlab("Gibbs free energy change of wild-type protein folding (in kcal/mol)") +
  ylab("Protein molecular weight (in Daltons)")




# Boxplot to gene groups - WORKS - keep as reference for now...


my_comparisons <- list( c("intracellular", "membrane"), c("intracellular", "secreted"), c("secreted", "membrane") )
ggboxplot(resGeneColl, x = "protType", y = "wtDG", color = "protType")+ 
  scale_color_manual(values=c("purple", "blue", "black")) +
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 6000)     # Add global p-value




# some quick checks, replace later
a <- subset(results, classificationVKGL == "LP")
b <- subset(results, classificationVKGL == "LB")
median(a$total.energy)
median(b$total.energy)
c <- subset(results, classificationVKGL == "LP" | classificationVKGL == "LB"| classificationVKGL == "VUS"  )

# by clsf
theme_set(theme_classic()) #or theme_bw()
ggplot(c, aes(x=total.energy, color=classificationVKGL, fill=classificationVKGL)) +
  theme_bw() +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=T) +
  scale_color_manual(values=c("black", "black", "black")) +
  scale_fill_manual(values=c("green", "red", "grey"))

# by chap
ggplot(c, aes(x=total.energy, color=chaperoned, fill=chaperoned)) +
  theme_bw() +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=T) +
  scale_color_manual(values=c("black", "black", "black")) +
  scale_fill_manual(values=c("yellow", "blue"))

youdenIndex = 5
tp <- sum(c[c$classificationVKGL=="LP",'total.energy'] >= youdenIndex)
fp <- sum(c[c$classificationVKGL=="LB",'total.energy'] >= youdenIndex)
tn <- sum(c[c$classificationVKGL=="LB",'total.energy'] < youdenIndex)
fn <- sum(c[c$classificationVKGL=="LP",'total.energy'] < youdenIndex)
ppv <- 100 *tp/(tp+fp)
npv <- 100 *tn/(tn+fn)
ppv
npv

ggplot(c, aes(x=protType, y=total.energy, fill=classificationVKGL)) +
  geom_violin()

d <- subset(results, classificationVKGL == "LP" | classificationVKGL == "LB")
ggplot(d %>% arrange(match(classificationVKGL, c("LB", "LP"))), aes(y=total.energy, x=wtDG, color=classificationVKGL)) +
  theme_classic() +
  geom_point(size=1) +
  scale_color_manual(values=c("green", "red", "grey")) +
  scale_x_continuous(labels = label_comma()) +
  xlab("Gibbs free energy change of wild-type protein folding (ΔG, in kcal/mol)") +
  ylab(paste("Difference in Gibbs free energy change between wild-type and variant protein (ΔΔG, in kcal/mol) ", sep=" "))

ggplot(d %>% arrange(match(classificationVKGL, c("LB", "LP"))), aes(y=total.energy, x=mwDa, color=classificationVKGL)) +
  theme_classic() +
  geom_point(size=1) +
  scale_color_manual(values=c("green", "red", "grey")) +
  scale_x_continuous(labels = label_comma()) +
  xlab("Protein molecular weight (in Daltons)") +
  ylab(paste("Difference in Gibbs free energy change between wild-type and variant protein (ΔΔG, in kcal/mol)", sep=" "))


