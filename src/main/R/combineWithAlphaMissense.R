library(seqminer)
library(toprdata)
library(ggplot2)
library(stringr)
library(scales)
library(cutpointr)
library(plyr)
library(dplyr)

rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability"
alphaMissenseLoc <- "/Applications/AlphaFold2/AlphaMissense_hg38.tsv.gz"
freeze1 <- paste(rootDir, "data", "freeze1.csv", sep="/")
results <- read.csv(freeze1)

uniqGenes <- unique(results$gene)
resultsWithAM <- data.frame()
for(i in 1:length(uniqGenes)){
  geneName <- uniqGenes[i]
  cat(paste("Working on ",geneName," (",i," of ",length(uniqGenes),")\n", sep=""))
  geneCoords <- subset(ENSGENES, gene_symbol==geneName)
  if(nrow(geneCoords)==0){
    cat(paste0("No coords for ", geneName, ", skipping..."))
    next
  }
  geneChr <- gsub("chr","", geneCoords$chrom)
  geneTabix <- paste(geneCoords$chrom, paste(geneCoords$gene_start, geneCoords$gene_end, sep="-"), sep=":")
  geneAlphaMissenseData <- tabix.read(alphaMissenseLoc, geneTabix)
  alphaMissense <- read.table(text=geneAlphaMissenseData, sep = '\t', header = FALSE, fileEncoding = "UTF-16LE", col.names = c("CHROM", "POS", "REF", "ALT", "genome", "uniprot_id", "transcript_id", "protein_variant", "am_pathogenicity", "am_class"))
  vkglGeneWithAlph <- merge(x = results, y = alphaMissense, by.x = c("dna_variant_pos", "dna_variant_ref", "dna_variant_alt"), by.y = c("POS","REF", "ALT"))
  if(!all(vkglGeneWithAlph$uniprot_id == vkglGeneWithAlph$UniProtID)){
    cat(paste0("UniProt ID mismatch for ", geneName, ", skipping..."))
    next
  }
  AAwithoutAchain <- paste0(substr(vkglGeneWithAlph$delta_aaSeq, 1, 1), substr(vkglGeneWithAlph$delta_aaSeq, 3, nchar(vkglGeneWithAlph$delta_aaSeq)))
  if(!all(AAwithoutAchain == vkglGeneWithAlph$protein_variant)){
    cat(paste0("AA change mismatch for ", geneName, ", skipping..."))
    next
  }
  cat(paste("Adding ",nrow(vkglGeneWithAlph)," variants\n", sep=""))
  resultsWithAM <- rbind(resultsWithAM, vkglGeneWithAlph)
}
# Persist these mappings to make it faster next time
resultsWithAMFile <- paste(rootDir, "data", "resultsWithAM.csv", sep="/")
#write.csv(resultsWithAM, resultsWithAMFile, row.names = FALSE, quote = FALSE)
resultsWithAM <- read.csv(resultsWithAMFile)


####################################################
# Determine optimal threshold using Youden's Index #
####################################################
# Read existing results
resultsForY <- subset(resultsWithAM, ann_classificationVKGL=="LB" | ann_classificationVKGL == "LP")
opt_cut <- cutpointr(resultsForY, am_pathogenicity, ann_classificationVKGL, direction = ">=", pos_class = "LP", neg_class = "LB", method = maximize_metric, metric = youden)
youdenIndex <- opt_cut$optimal_cutpoint
tp <- sum(resultsForY[resultsForY$ann_classificationVKGL=="LP",'am_pathogenicity'] >= youdenIndex)
fp <- sum(resultsForY[resultsForY$ann_classificationVKGL=="LB",'am_pathogenicity'] >= youdenIndex)
tn <- sum(resultsForY[resultsForY$ann_classificationVKGL=="LB",'am_pathogenicity'] < youdenIndex)
fn <- sum(resultsForY[resultsForY$ann_classificationVKGL=="LP",'am_pathogenicity'] < youdenIndex)
ppv <- 100 *tp/(tp+fp)
npv <- 100 *tn/(tn+fn)
sens <- opt_cut$sensitivity*100
spec <- opt_cut$specificity*100
cat(paste("The optimal AlphaMissense threshold based on ",dim(results)[1]," VKGL variant classifications is ",round(youdenIndex,5)," with PPV ",round(ppv),"%, NPV ",round(npv),"%, sensitivity ",round(sens),"% and specificity ",round(spec),"%.\n",sep=""))


#################
# Other outputs #
#################
# Confusion matrix between AlphaMissense label and VKGL label
table(resultsForY$am_class, resultsForY$ann_classificationVKGL)
# Density plot of AlphaMissense predictions vs VKGL label
ggplot(resultsForY, aes(am_pathogenicity, colour = ann_classificationVKGL, fill = ann_classificationVKGL)) +
  theme_bw() + theme(panel.grid = element_blank(), axis.title.x=element_text(size=10)) +
  geom_density(alpha = 0.1, adjust = 0.1) +
  scale_fill_manual(values=c("green", "red")) +
  scale_color_manual(values=c("green", "red")) +
  geom_vline(xintercept = youdenIndex) +
  ggtitle(paste("Optimal AlphaMissense threshold based on\nVKGL variant classifications: ",round(youdenIndex,5),sep="")) +
  xlab("AlphaMissense pathogenicity prediction") +
  ylab("Density") +
  guides(fill=guide_legend(title="VKGL\nvariant\nclassifi-\ncation")) +
  guides(color=guide_legend(title="VKGL\nvariant\nclassifi-\ncation"))
#ggsave(paste("alphamissense-vkgl-threshold-density.png",sep=""), width=5, height=5)


###############################
# RF with AlphaMissense added #
###############################

# AUC and plot for only AlphaMissense pathogenicity scores
rf.pred = prediction(resultsForY$am_pathogenicity, resultsForY$ann_classificationVKGL)
rf.perf = performance(rf.pred, "tpr", "fpr")
auc <- performance(rf.pred,"auc")
auc <- unlist(slot(auc, "y.values"))
plot(rf.perf,main=paste0("Variant classification on AlphaMissense, RF ROC curve (AUC ",round(auc,2),")"),col=2,lwd=2)

# Combine peptide, FoldX and AlphaMissense in Random Forest model
set.seed(222)
resultsForY$ann_classificationVKGL <- as.factor(resultsForY$ann_classificationVKGL)
draw <- sample(c(TRUE, FALSE), nrow(resultsForY), replace=TRUE, prob=c(0.8, 0.2))
train <- resultsForY[draw, ]
test <- resultsForY[!draw, ]

rf <-randomForest(ann_classificationVKGL~., data=train, mtry=7, ntree=1000, keep.forest=TRUE, importance=TRUE, xtest=subset(test, select=-ann_classificationVKGL))
varImpPlot(rf)
rf.pr = predict(rf,type="prob",newdata=test)[,2]
rf.pred = prediction(rf.pr, test$ann_classificationVKGL)
rf.perf = performance(rf.pred, "tpr", "fpr")

auc <- performance(rf.pred,"auc")
auc <- unlist(slot(auc, "y.values"))
auc
plot(rf.perf,main=paste0("Variant classification on peptide and folding properties\nplus AlphaMissense, RF ROC curve (AUC ",round(auc,2),")"),col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")


