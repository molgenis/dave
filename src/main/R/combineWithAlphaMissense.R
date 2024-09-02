library(seqminer)
library(toprdata)
library(ggplot2)
library(stringr)
library(scales)
library(cutpointr)
library(plyr)
library(dplyr)
library(crunch)    # Compress results

rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability"
alphaMissenseLoc <- "/Applications/AlphaFold2/AlphaMissense_hg38.tsv.gz"
freeze2 <- paste(rootDir, "data", "freeze2.csv.gz", sep="/")
results <- read.csv(freeze2)
results$AAwithoutAchain <- paste0(substr(results$delta_aaSeq, 1, 1), substr(results$delta_aaSeq, 3, nchar(results$delta_aaSeq)))

uniqGenes <- unique(results$gene)
resultsWithAM <- data.frame()
for(i in 1:length(uniqGenes)){
  geneName <- uniqGenes[i]
  variantsForGene <- subset(results, gene == geneName)
  chrom <- unique(variantsForGene$dna_variant_chrom)
  if(length(chrom) > 1)
  {
    stop(paste0("Multiple chroms for ", geneName, ", stopping"))
  }
  chrom <- paste0("chr", chrom)
  gene_start <- min(variantsForGene$dna_variant_pos)
  gene_end <- max(variantsForGene$dna_variant_pos)
  
  cat(paste("Working on ",geneName," (",i," of ",length(uniqGenes),")\n", sep=""))
  geneTabix <- paste(chrom, paste(gene_start, gene_end, sep="-"), sep=":")
  geneAlphaMissenseData <- tabix.read(alphaMissenseLoc, geneTabix)
  if(length(geneAlphaMissenseData)==0)
  {
    cat(paste0("No AlphaMissense data for ", geneName, ", skipping"))
    next
  }
  alphaMissense <- read.table(text=geneAlphaMissenseData, sep = '\t', header = FALSE, fileEncoding = "UTF-16LE", colClasses = c("character", "numeric", "character", "character", "character", "character", "character", "character", "numeric", "character"), col.names = c("CHROM", "POS", "REF", "ALT", "genome", "uniprot_id", "transcript_id", "protein_variant", "am_pathogenicity", "am_class"))
  alphaMissense$CHROM <- gsub("chr", "", alphaMissense$CHROM)
  # Sometimes a merge on seemly good data results in 0 overlap, e.g. for ALPK3
  # On close inspection, the AA change is different, probably due to trancript subversions (ENST00000258888 vs. AM: ENST00000258888.5)
  # Also, we could match on transcript here as well, but that would require a 'contains' operation (see later)
  vkglGeneWithAlph <- merge(x = variantsForGene, y = alphaMissense, by.x = c("dna_variant_chrom", "dna_variant_pos", "dna_variant_ref", "dna_variant_alt","AAwithoutAchain","UniProtID"), by.y = c("CHROM","POS","REF", "ALT","protein_variant","uniprot_id"))
  if(dim(vkglGeneWithAlph)[1]==0)
  {
    cat(paste0("No merged data for ", geneName, ", skipping"))
    next
  }
  vkglGeneWithAlph$key <- paste0(vkglGeneWithAlph$dna_variant_chrom,"_",vkglGeneWithAlph$dna_variant_pos,"_",vkglGeneWithAlph$dna_variant_ref,"_",vkglGeneWithAlph$dna_variant_alt,"_",vkglGeneWithAlph$AAwithoutAchain,"_",vkglGeneWithAlph$UniProtID)
  #cat(paste("Adding ",nrow(vkglGeneWithAlph)," variants\n", sep=""))
  resultsWithAM <- rbind(resultsWithAM, vkglGeneWithAlph)
}
# The merged results contain some overlapping genes, cannot use directly
results$ann_am_pathogenicity <- NA
for(i in 1:nrow(results)){
  #i <- 3523
  row <- results[i,]
  resultKey <- paste0(row$dna_variant_chrom,"_",row$dna_variant_pos,"_",row$dna_variant_ref,"_",row$dna_variant_alt,"_",row$AAwithoutAchain,"_",row$UniProtID)
  fromAM <- resultsWithAM[resultsWithAM$key==resultKey,]
  if(dim(fromAM)[1]==0){
    #cat(paste0("No data for ", i, ", skipping\n"))
  }else if(dim(fromAM)[1]==1){
    results[i,"ann_am_pathogenicity"] <- fromAM$am_pathogenicity
  }else{
    # Sometimes the exact same variant was present for multiple transcript, e.g. 19_53982805_G_A_A412T_Q8WXS5 for ENST00000270458.2 and ENST00000638874.1.
    # However, AM pathogenicity score is the same, in this case 0.1358.
    #cat(paste0("Multiple data for ", i, ", skipping\n"))
  }
}
# Remove tmp column used for merging
results$AAwithoutAchain <- NULL

# Persist these mappings to make it faster next time
resultsWithAMFile <- paste(rootDir, "data", "freeze2-AMmerge.csv.gz", sep="/")
#write.csv.gz(results, resultsWithAMFile, row.names = FALSE, quote = FALSE)
resultsWithAM <- read.csv(resultsWithAMFile)


####################################################
# Determine optimal threshold using Youden's Index #
####################################################
# Read existing results
resultsForY <- subset(resultsWithAM, (ann_classificationVKGL=="LB" | ann_classificationVKGL == "LP") ) #  &  (ann_proteinLocalization == "membrane" | ann_proteinLocalization == "intracellular") / ann_proteinLocalization == "secreted"
resultsForY <- resultsForY[!is.na(resultsForY$ann_am_pathogenicity), ]
resultsForY$ann_classificationVKGL <- as.factor(resultsForY$ann_classificationVKGL)
opt_cut <- cutpointr(resultsForY, ann_am_pathogenicity, ann_classificationVKGL, direction = ">=", pos_class = "LP", neg_class = "LB", method = maximize_metric, metric = youden)
opt_cut$AUC
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


