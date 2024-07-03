library(ggplot2)   # for plotting
library(ggpubr)    # for boxplots with pvalues
library(dplyr)     # to remove duplicate rows
library(scales)    # for big values with commas in plots

rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability"
freeze1 <- paste(rootDir, "data", "freeze1.csv", sep="/")
results <- read.csv(freeze1)


################################################################################
# Collapse variant-level results into gene-level and make scatterplot
#######################################################################################
resGeneColl <- data.frame(gene=results$gene, mwDa=results$mwDa, wtDG=results$wtDG, transcript=results$transcript, uniprot=results$uniprot, protType=results$protType, chaperoned=results$chaperoned)
resGeneColl <- resGeneColl %>% distinct()
resGeneColl$energyPerKDa <- resGeneColl$wtDG/(resGeneColl$mwDa/1000)
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
ggboxplot(resGeneColl, x = "protType", y = "mwDa", color = "protType")+ 
  scale_color_manual(values=c("purple", "blue", "black")) +
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 6000)     # Add global p-value

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


