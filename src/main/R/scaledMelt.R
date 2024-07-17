library(scales)
library(reshape)
library(reshape2)


rootDir <- "C:/Users/tk_20/git/vkgl-secretome-protein-stability"
freeze1 <- paste(rootDir, "data", "freeze1.csv", sep="/")
results <- read.csv(freeze1)
rownames(results) <- paste0(results$gene, ":", results$delta_aaSeq)

nums <- unlist(lapply(results, is.numeric), use.names = FALSE)
results <- results[ , nums]

results <- results[, c("mutant_bomanIndex", "mutant_sloop_entropy", "mutant_instabilityIndex")]

results_scaled <- scale(results)
results_sc_melt <- reshape2::melt(results_scaled,  na.rm = FALSE, value.name = "ScaledValue", varnames = c("GeneMut", "Property"), id=c("mutant_bomanIndex","mutant_sloop_entropy", "mutant_instabilityIndex"))
head(results_sc_melt)


