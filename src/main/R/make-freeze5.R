######################
# Load used packages #
######################
library(crunch)     # To compress results

#######################
# Adjustable settings #
#######################
vkglProtVarFileName <- "VKGL_apr2024_protForFolding.tsv" # loading processed VKGL protein variants from: /data/genes/{gene}/{vkglProtVarFile}
# Example for Unix
rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability" # root directory that contains README.md, data/, img/, out/, src/, etc.
# Example for Windows
#rootDir <- "D:/github/vkgl-secretome-protein-stability"

# Assuming data was produced by main.R and enriched by CombineWithAlhaMissense.R
# Load the freeze4 data as starting point for creating freeze5
freeze4loc <- paste(rootDir, "data", "freeze4.csv.gz", sep="/")
freeze4 <- read.csv(freeze4loc)
rownames(freeze4) <- paste0(freeze4$gene, "/", freeze4$UniProtID, ":", freeze4$delta_aaSeq)
dataGenesDir <- paste(rootDir, "data", "genes", sep="/")
source(paste(rootDir, "src", "main", "R", "helper-functions.R", sep="/"))
colnames(freeze4)

# Where applicable, use delta's only (= difference between WT and mutant)
subFreeze4 <- freeze4[, c(
                      # Basic information: gene annotations, genomic position, DNA and protein mutation (n = 9)
                      "gene", "TranscriptID", "UniProtID", "dna_variant_chrom", "dna_variant_pos", "dna_variant_ref", "dna_variant_alt", "dna_variant_assembly", "delta_aaSeq",
                      # VKGL classification label (n = 1)
                      "ann_classificationVKGL",
                      # Protein annotations (n = 2)
                      "ann_proteinIschaperoned", "ann_proteinLocalization",
                      # FoldX energy terms (n = 22)
                      "delta_total.energy", "delta_Backbone.Hbond", "delta_Sidechain.Hbond", "delta_Van.der.Waals", "delta_Electrostatics", "delta_Solvation.Polar", "delta_Solvation.Hydrophobic", "delta_Van.der.Waals.clashes", "delta_entropy.sidechain", "delta_entropy.mainchain", "delta_sloop_entropy", "delta_mloop_entropy", "delta_cis_bond", "delta_torsional.clash", "delta_backbone.clash", "delta_helix.dipole", "delta_water.bridge", "delta_disulfide", "delta_electrostatic.kon", "delta_partial.covalent.bonds", "delta_energy.Ionisation", "delta_Entropy.Complex", 
                      # Other relevant protein properties (n = 9)
                      "delta_aliphaticIndex", "delta_bomanIndex", "delta_charge", "delta_hydrophobicMoment", "delta_hydrophobicity", "delta_instabilityIndex", "delta_molWeight", "delta_massOverCharge", "delta_isoElecPoint",
                      # AlphaMissense pathogenicity prediction, for reference (n = 1)
                      "ann_am_pathogenicity")
                      # TOTAL: n = 44
                   ]
colnames(subFreeze4)
freeze4 <- NULL

# Impute AlphaMissense in source data with mean
am_mean <- mean(subFreeze4[!is.na(subFreeze4$ann_am_pathogenicity),]$ann_am_pathogenicity)
subFreeze4["ann_am_pathogenicity"][is.na(subFreeze4["ann_am_pathogenicity"])] <- am_mean

# Genes present in freeze 4, meaning they have at least 1 succesfully analyzed mutation
succesfulGenes <- unique(subFreeze4$gene)

# TEMPORARY: get imputed values for GeoNet per classification label or globally (!)
#gnImp <- get_GeoNet_bellcurve_values(dataGenesDir, vkglProtVarFileName, succesfulGenes, FALSE)

# Iterate over genes and gather functional features
funcFeat <- data.frame()
for(i in seq_along(succesfulGenes))
{
  # i <- 1 # DEBUG/DEV
  geneName <- succesfulGenes[i]
  cat(paste("Working on gene:", geneName, "\n", sep=" "))
  specificGeneDir <- paste(dataGenesDir, geneName, sep="/")
  
  glmScore_WT <- read.csv(paste(specificGeneDir, "DNA_RNA_interaction_terms.csv", sep="/"))
  p2rank_WT <- read.csv(paste(specificGeneDir, list.files(specificGeneDir, pattern="*_v4_Repair.pdb_predictions.csv"), sep="/"))
  
  imputeGN <- FALSE
  if(length(list.files(specificGeneDir, pattern="CombinedGeoNetPredictions.csv")) == 0){
    imputeGN <- TRUE
  }else{
    geoNet_WT <- read.csv(paste(specificGeneDir, "CombinedGeoNetPredictions.csv", sep="/"))
  }
  
  variants <- read.table(file=paste(specificGeneDir, vkglProtVarFileName, sep="/"), sep = '\t', header = TRUE, colClasses = c("character", "character", "numeric", "character", "character", "character", "character", "numeric", "character"))
  foldingResultsDir <- paste(specificGeneDir, "folding-results", sep="/")
  
  for(j in 1:nrow(variants))
  {
    # j <- 1 # DEBUG/DEV
    mutation <- variants$ProtChange[j]
    cat(paste("Getting results from ", mutation, " (gene ",geneName,", mutation ", j, " of ", nrow(variants), ")\n", sep=""))
    mutationDir <- paste(foldingResultsDir, mutation, sep="/")
    
    if(length(list.files(mutationDir, pattern="exception.txt")) > 0)
    {
      cat("  exception found, skipping......\n")
      next
    }
    
    glmScore_Mu <- read.csv(paste(mutationDir, "DNA_RNA_interaction_terms.csv", sep="/"))
    p2rank_Mu <- read.csv(paste(mutationDir, list.files(mutationDir, pattern="*_v4_Repair.pdb_predictions.csv"), sep="/"))
    
    glmFeat <- extract_glmscore_features(glmScore_WT, glmScore_Mu)
    p2rankFeat <- extract_p2rank_features(p2rank_WT, p2rank_Mu)
    
    if(!imputeGN && !length(list.files(mutationDir, pattern="CombinedGeoNetPredictions.csv")) == 0){
      geoNet_Mu <- read.csv(paste(mutationDir, "CombinedGeoNetPredictions.csv", sep="/"))
      geonetFeat <- extract_geonet_features(geoNet_WT, geoNet_Mu)
    }else{
      #geonetFeat <- impute_GeoNet_features(gnImp, variants$Classification[j])
      geonetFeat <- na_GeoNet_features()
    }
    funcFeat <- rbind(funcFeat, data.frame(p2rankFeat, glmFeat, geonetFeat))
  }
}

# Combine into freeze 5
freeze5raw <- cbind(subFreeze4, funcFeat)
dim(freeze5raw)

# Remove all columns with only 0 values
freeze5 <- freeze5raw[, colSums(freeze5raw != 0, na.rm = T) > 0]
dim(freeze5)

# Inspect removed columns
setdiff(names(freeze5raw), names(freeze5))

# Save as freeze 5
resultsFreeze5 <- paste(rootDir, "data", "freeze5-provisional-NAs.csv.gz", sep="/")
write.csv.gz(freeze5, resultsFreeze5, row.names = FALSE, quote = FALSE)

# Load back in
#frz5 <- read.csv(resultsFreeze5)
#rownames(frz5) <- paste0(frz5$gene, "/", frz5$UniProtID, ":", frz5$delta_aaSeq)
#head(frz5)
