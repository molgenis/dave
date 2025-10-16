library(apcluster)
library(dplyr)

rootDir <- "/Users/joeri/git/dave1"

# Variants, slow
freeze1 <- paste(rootDir, "data", "freeze1.csv", sep="/")
data <- read.csv(freeze1)
rownames(data) <- paste0(data$gene, ":", data$delta_aaSeq)
# Subset
data <- subset(data, ann_classificationVKGL == "LP") # just LP
# Remove various columns that are not so relevant
nums <- unlist(lapply(data, is.numeric), use.names = FALSE)
numData <- data[ , nums]
numData <- subset(numData, select = -c(dna_variant_pos, ann_mutant_energy_SD, ann_WT_energyPerKDa))
numData <- numData %>% select(-contains("WT_"))
# Run APclust
numDataNegDist <- negDistMat(numData, r=2)
apclust <- apcluster(numDataNegDist)
apclust
plot(apclust, numData[, c(1, 2, 3, 4, 5, 6)])

# Proteins, easy
peptidePropPerGeneFile <- paste(rootDir, "data", "peptidePropPerGene.csv", sep="/")
peptidePropPerGene <- read.csv(peptidePropPerGeneFile)
rownames(peptidePropPerGene) <- peptidePropPerGene$WT_gene
nums <- unlist(lapply(peptidePropPerGene, is.numeric), use.names = FALSE)
numDat <- peptidePropPerGene[ , nums]
numDatNegDist <- negDistMat(numDat, r=2)
apclust <- apcluster(numDatNegDist)
apclust
plot(apclust, numDat[, c(1, 2, 3, 4, 5, 6)])
