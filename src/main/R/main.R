######################
# Load used packages #
######################
library(R.utils)   # for 'gunzip', 'mkdirs'
library(Rpdb)      # to load PDB files
library(dplyr)     # to remove duplicate rows
library(Peptides)  # protein annotations
library(peptoolkit)# QSAR features
library(crunch)    # Compress results
library(protr)     # PseAAC/APseAAC features


#######################
# Adjustable settings #
#######################
# Example for Unix
rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability" # root directory that contains README.md, data/, img/, out/, src/, etc.
vkglProtVarFileName <- "VKGL_apr2024_protForFolding.tsv" # loading processed VKGL protein variants from: /data/genes/{gene}/{vkglProtVarFile}
foldxExec <- "/Applications/FoldX/5/foldx5MacStd/foldx_20241231" # exact path to the FoldX executable
alphaFoldLoc <- "/Applications/AlphaFold2/mane_overlap_v4.tar"
UP000005640_9606_HUMAN_v4Loc <- "/Applications/AlphaFold2/UP000005640_9606_HUMAN_v4.tar"
# Example for Windows
rootDir <- "D:/github/vkgl-secretome-protein-stability"
vkglProtVarFileName <- "VKGL_apr2024_protForFolding.tsv"
foldxExec <- "C:/\"Program Files\"/FoldX/foldx_20241231.exe"
alphaFoldLoc <- "D:/mane_overlap_v4.tar"
UP000005640_9606_HUMAN_v4Loc <- "D:/UP000005640_9606_HUMAN_v4.tar"
# Possible gene lists to work on
secr <- read.table(file=paste(rootDir, "data", "protein-atlas-secreted-geneIDs-mane-uniprot-withvariants.txt", sep="/"), sep = '\t',header = TRUE)
secr$protType <- "secreted"
intr <- read.table(file=paste(rootDir, "data", "protein-atlas-intracellular-geneIDs-mane-uniprot-random2000-withvariants.txt", sep="/"), sep = '\t',header = TRUE)
intr$protType <- "intracellular"
memb <- read.table(file=paste(rootDir, "data", "protein-atlas-membrane-geneIDs-mane-uniprot-random2000-withvariants.txt", sep="/"), sep = '\t',header = TRUE)
memb$protType <- "membrane"
allg <- rbind(secr, intr, memb)
# Selected gene list to work on
selectedGenes <- allg


##########################################
# Derived paths of directories and files #
##########################################
dataGenesDir <- paste(rootDir, "data", "genes", sep="/")
source(paste(rootDir, "src", "main", "R", "helper-functions.R", sep="/"))


#################################################
# Sort the mapping of HGNC symbol to UniProt ID #
#################################################
head(selectedGenes)
geneNames <- sort(selectedGenes$Gene.name)
head(geneNames)


######################################################
# Iterate over all gene directories, repair and fold #
######################################################
setwd(dataGenesDir)
for(i in seq_along(geneNames))
{
  geneName <- geneNames[i]
  cat(paste("Working on ", geneName, " (gene ", i, " of ", length(geneNames), ")\n", sep=""))
  specificGeneDir <- paste(dataGenesDir, geneName, sep="/")
  setwd(specificGeneDir)
  if(length(list.files(pattern="exception.txt")) > 0)
  {
    cat("  gene already tried before but failed, continuing...\n")
    next
  }
  uniProtID <- selectedGenes$UniProtKB.Swiss.Prot.ID[selectedGenes$Gene.name==geneName]
  # first cleanup repair fxout if present
  file.remove(list.files(specificGeneDir, pattern="*.fxout"))
  file.remove(list.files(specificGeneDir, pattern="molecules", include.dirs=TRUE))
  file.remove(list.files(pattern="rotabase.txt"))
  # Find in TAR file, extract into working dir, gunzip and repair (may take a while)
  if(length(list.files(pattern="*_Repair.pdb")) == 0){
    alphaFoldAll <- untar(alphaFoldLoc, list = TRUE)
    alphaFoldPDBs <- grep(".pdb.gz",alphaFoldAll, value=TRUE)
    PDBForGeneGz <- grep(paste(uniProtID,collapse="|"), alphaFoldPDBs, value=TRUE)
    if(identical(PDBForGeneGz, character(0))){
      cat("  gene not found in MANE AlphaFold, trying UP000005640_9606_HUMAN_v4...\n")
      UP000005640_9606_HUMAN_v4_All <- untar(UP000005640_9606_HUMAN_v4Loc, list = TRUE)
      UP000005640_9606_HUMAN_v4_PDBs <- grep(".pdb.gz", UP000005640_9606_HUMAN_v4_All, value=TRUE)
      PDBForGeneGz <- grep(paste(uniProtID,collapse="|"), UP000005640_9606_HUMAN_v4_PDBs, value=TRUE)
      if(identical(PDBForGeneGz, character(0))){
        write(paste("No PDB found for uniprotID", uniProtID), file = "exception.txt")
        next
      }else{
        untar(UP000005640_9606_HUMAN_v4Loc, files = PDBForGeneGz)
      }
    }else{
      untar(alphaFoldLoc, files = PDBForGeneGz)
    }
    if(length(PDBForGeneGz) > 1){
      write(paste("Multiple PDB and/or fragments found for uniprotID", uniProtID, ":", PDBForGeneGz), file = "exception.txt")
      file.remove(list.files(pattern="*.pdb.gz"))
      next
    }
    PDBForGene <- gunzip(PDBForGeneGz, overwrite=TRUE)[[1]]
    system(paste(foldxExec, " --command=RepairPDB --pdb=",PDBForGene,sep=""), intern = TRUE)
  } else{
    # Verify that uniprot ID matched currently repaired protein
    if(!grepl(uniProtID, list.files(pattern="*_Repair.pdb"), fixed=TRUE)){
      stop(paste0("Repair file based on ", uniProtID, " expected but found ", list.files(pattern="*_Repair.pdb")))
    }
    # Cleanup original PDB if a repaired PDB is present
    file.remove(list.files(pattern="*model_v4.pdb"))
  }
  repPDB_fileName <- list.files(pattern="*_Repair.pdb")
  repPDB_fullLoc <- paste(specificGeneDir, repPDB_fileName, sep="/")
  repPDB_fullLoc
  
  # Read the VKGL mutations for this protein
  vkgl <- read.table(file=paste(specificGeneDir, vkglProtVarFileName, sep="/"), sep = '\t', header = TRUE)
  mutations <- vkgl$ProtChange
  
  # Calculate change in protein stability for each mutation
  foldingResultsDir <- paste(specificGeneDir, "folding-results", sep="/")
  mkdirs(foldingResultsDir)
  for(j in seq_along(mutations))
  {
    cat(paste("Working on ", mutations[j], " (gene ",geneName,", mutation ", j, " of ", length(mutations), ")\n", sep=""))
    mutationDir <- paste(foldingResultsDir, mutations[j], sep="/")
    mkdirs(mutationDir)
    setwd(mutationDir)
    file.remove(list.files(pattern="molecules", include.dirs=TRUE))
    file.remove(list.files(pattern="rotabase.txt"))
    file.remove(list.files(pattern="*.pdb"))
    file.remove(list.files(pattern="*individual_list.txt"))
    file.remove(list.files(pattern="PdbList_*"))
    file.remove(list.files(pattern="Dif_*"))
    file.remove(list.files(pattern="Raw_*"))
    if(length(list.files(mutationDir, pattern="*.fxout")) > 0)
    {
      cat("  already folded, skipping...\n")
      next
    }
    if(length(list.files(mutationDir, pattern="exception.txt")) > 0)
    {
      cat("  already tried before but failed, skipping......\n")
      next
    }
    file.copy(from = repPDB_fullLoc, to = mutationDir)
    write(paste(vkgl[j, "ProtChange"], ";", sep=""), file = "individual_list.txt")
    state <- system(paste(foldxExec, " --command=BuildModel --mutant-file=individual_list.txt --numberOfRuns=5 --pdb=", repPDB_fileName, sep=""), intern = TRUE)
    if(any(grepl("Specified residue not found", state)))
      {
        write(state, file = "exception.txt")
        #break
    }
     file.remove(repPDB_fileName)
     cat("...done!\n")
  }
}


###############################
# Add wild type info for gene #
###############################
setwd(dataGenesDir)
for(i in seq_along(geneNames))
{
  geneName <- geneNames[i]
  specificGeneDir <- paste(dataGenesDir, geneName, sep="/")
  setwd(specificGeneDir)
  pdbFile <- list.files(pattern="*_Repair.pdb")
  if(length(list.files(pattern="wt.txt")) == 0 & length(pdbFile) > 0){
    vkgl <- read.table(file=paste(specificGeneDir, vkglProtVarFileName, sep="/"), sep = '\t', header = TRUE)
    mutationsToCheck <- vkgl$ProtChange
    successfulMut <- NULL
    for(mutation in mutationsToCheck){
      mutationResultDir <- paste(specificGeneDir, "folding-results", mutation, sep="/")
      avgEnergy <- list.files(mutationResultDir, pattern="Average")
      if(length(avgEnergy) > 0){
        successfulMut <- mutation
        break
      }
    }
    if(is.null(successfulMut))
    {
      cat(paste("No average energy file for gene, skipping", geneName, "\n", sep=" "))
      next
    }
      cat(paste("Gene", geneName, "has at least 1 successful mutation:", successfulMut,"so doing quick refold to get WT DG...\n", sep=" "))
      tmpDir <- paste(specificGeneDir, "tmp", sep="/")
      mkdirs(tmpDir)
      setwd(tmpDir)
      file.copy(from = paste(specificGeneDir, pdbFile, sep="/"), to = tmpDir)
      write(paste(successfulMut, ";", sep=""), file = "individual_list.txt")
      state <- system(paste(foldxExec, " --command=BuildModel --mutant-file=individual_list.txt --numberOfRuns=1 --pdb=", pdbFile, sep=""), intern = TRUE)
      rawEnergy <- list.files(pattern="Raw")
      rawResult <- read.table(file = rawEnergy, header = TRUE, skip = 8, sep="\t")
      WTrow <- rawResult[2,]
      if(grepl("WT_AF-", WTrow$Pdb, fixed=TRUE)){
        write.csv(WTrow, file = paste(specificGeneDir, "wt.txt", sep="/"), row.names = FALSE, quote = FALSE)
        file.remove(list.files(pattern="*"), include.dirs=TRUE)
        setwd(specificGeneDir)
        file.remove(list.files(pattern="tmp"), include.dirs=TRUE)
        cat(paste("Done! Added wild type info for gene", geneName, "\n", sep=" "))
      }else{
        stop(paste0("Something went wrong with gene ", geneName, ", mutation ", + successfulMut))
      }
  }else{
    cat(paste("No PDB file or wild type info already present for gene", geneName, "\n", sep=" "))
  }
}


##################
# Gather results #
##################
setwd(dataGenesDir)
results <- data.frame()
geneLvlData <- data.frame()
for(i in seq_along(geneNames))
{
  geneName <- geneNames[i]
  geneInfo <- selectedGenes[selectedGenes$Gene.name==geneName, ]
  cat(paste("Loading data for ", geneName, " (gene ", i, " of ", length(geneNames), ")\n", sep=""))
  specificGeneDir <- paste(dataGenesDir, geneName, sep="/")
  setwd(specificGeneDir)
  wtInfoFile <- list.files(pattern="wt.txt")
  if(length(wtInfoFile) > 0){
    wtInfo <- read.csv(file = wtInfoFile, header = TRUE)
    geneInfo <- cbind(geneInfo, wtInfo[1,])
    geneLvlData <- rbind(geneLvlData, geneInfo)
  }
  variants <- read.table(file=paste(specificGeneDir, vkglProtVarFileName, sep="/"), sep = '\t', header = TRUE, colClasses = c("character", "character", "numeric", "character", "character", "character", "character", "numeric", "character"))
  foldingResultsDir <- paste(specificGeneDir, "folding-results", sep="/")
  for(j in 1:nrow(variants))
  {
    mutationDir <- paste(foldingResultsDir, variants$ProtChange[j], sep="/")
    if(!dir.exists(mutationDir))
    {
      next
    }
    if(length(list.files(mutationDir, pattern="*.fxout")) == 0){
      next
    }
    avgDiff <- list.files(mutationDir, pattern="Average")
    result <- read.table(file = paste(mutationDir, avgDiff, sep="/"), header = TRUE, skip = 8, sep="\t")
    result$assembly <- variants[j, "Assembly"]
    result$chrom <- variants[j, "Chrom"]
    result$pos <- variants[j, "Pos"]
    result$ref <- variants[j, "Ref"]
    result$alt <- variants[j, "Alt"]
    result$gene <- geneName
    result$protChange <- variants[j, "ProtChange"]
    result$classificationVKGL <- variants[j, "Classification"]
    results <- rbind(results, result)
  }
}


#####################################################################################
# Get and store all AA sequences for each gene with 1+ succesfully folded mutations #
#####################################################################################
setwd(dataGenesDir)
succesfulGenes <- unique(results$gene)
peptidePropPerGene <- data.frame(gene=character(), aliphaticIndex=numeric(), bomanIndex=numeric(), charge=numeric(), hydrophobicMoment=numeric(), hydrophobicity=numeric(), instabilityIndex=numeric(), molWeight=numeric(), massOverCharge=numeric(), isoElecPoint=numeric(),  aaSeq=character())
for(i in seq_along(succesfulGenes))
{
  geneName <- succesfulGenes[i]
  cat(paste("retrieving AA seq for gene: ", geneName, "\n", sep=" "))
  specificGeneDir <- paste(dataGenesDir, geneName, sep="/")
  setwd(specificGeneDir)
  pdbFile <- list.files(pattern="*_Repair.pdb")
  if(length(pdbFile) > 0){
    x <- read.pdb(pdbFile)
    collapsePDB <- data.frame(resid=x$atoms$resid, resname=x$atoms$resname)
    collapsePDB <- collapsePDB %>% distinct()
    aaSeq <- aa3to1(collapsePDB$resname)
    peptidePropPerGene <- rbind(peptidePropPerGene, data.frame(gene = geneName,
                                             aliphaticIndex = aIndex(aaSeq),
                                             bomanIndex = boman(aaSeq),
                                             charge = charge(aaSeq),
                                             hydrophobicMoment = hmoment(aaSeq),
                                             hydrophobicity = hydrophobicity(aaSeq),
                                             instabilityIndex = instaIndex(aaSeq),
                                             molWeight = mw(aaSeq),
                                             massOverCharge = mz(aaSeq),
                                             isoElecPoint = pI(aaSeq),
                                             aaSeq = aaSeq
                                             ))
    cat(paste0("calculations done\n"))
  }else{
    cat(paste("No PDB file for gene", geneName, "\n", sep=" "))
  }
}
# Merge with chaperone data
chapInteractingGenesLoc <- paste(rootDir, "data", "chaperones", "interacting-with-chaperones-genenames-merged-with-uniprotmapped.txt", sep="/")
chapInteractingGenes <- read.table(file=chapInteractingGenesLoc, sep = '\t',header = TRUE)
peptidePropPerGene$chaperoned <- as.factor(peptidePropPerGene$gene %in% chapInteractingGenes$Gene.name)
# Add wild-type folding data including transcript/uniprot/localization 
peptidePropPerGene <- merge(peptidePropPerGene, geneLvlData, by.x = "gene", by.y = "Gene.name")
peptidePropPerGene$energyPerKDa <- peptidePropPerGene$total.energy/(peptidePropPerGene$molWeight/1000)
# Rename columns to WT ("Wild Type") for clarity
names(peptidePropPerGene)[names(peptidePropPerGene) == 'protType'] <- 'proteinLocalization'
colnames(peptidePropPerGene) <- paste('WT', colnames(peptidePropPerGene), sep = '_')
# Persist these results for quick loading later
peptidePropPerGeneFile <- paste(rootDir, "data", "peptidePropPerGene.csv", sep="/")
#write.csv(peptidePropPerGene, peptidePropPerGeneFile, row.names = FALSE, quote = FALSE)
peptidePropPerGene <- read.csv(peptidePropPerGeneFile)


################################################################################
# For all succesfully folded mutations, mutate AA-seq and calculate properties #
################################################################################
setwd(dataGenesDir)
succesfulGenes <- unique(peptidePropPerGene$WT_gene)
results[,"mutant_aliphaticIndex"] <- NA
results[,"mutant_bomanIndex"] <- NA
results[,"mutant_charge"] <- NA
results[,"mutant_hydrophobicMoment"] <- NA
results[,"mutant_hydrophobicity"] <- NA
results[,"mutant_instabilityIndex"] <- NA
results[,"mutant_molWeight"] <- NA
results[,"mutant_massOverCharge"] <- NA
results[,"mutant_isoElecPoint"] <- NA
results[,"delta_aliphaticIndex"] <- NA
results[,"delta_bomanIndex"] <- NA
results[,"delta_charge"] <- NA
results[,"delta_hydrophobicMoment"] <- NA
results[,"delta_hydrophobicity"] <- NA
results[,"delta_instabilityIndex"] <- NA
results[,"delta_molWeight"] <- NA
results[,"delta_massOverCharge"] <- NA
results[,"delta_isoElecPoint"] <- NA
results[,"mutant_aaSeq"] <- NA
for(i in 1:nrow(results))
{
  variant <- results[i,]
  cat(paste0("Calculating variant peptide properties for ", variant$gene, ":", variant$protChange,"\n"))
  genePeptideProp <- peptidePropPerGene[peptidePropPerGene$WT_gene==variant$gene, ]
  aaToReplace <- substr(variant$protChange, 1, 1)
  aaReplacement <- substr(variant$protChange, nchar(variant$protChange), nchar(variant$protChange))
  aaReplacePos <- as.numeric(substr(variant$protChange, 3, (nchar(variant$protChange)-1)))
  aaSeqToMutate <- genePeptideProp$WT_aaSeq
  originalAAatPos <- substr(aaSeqToMutate, aaReplacePos, aaReplacePos)
  if(originalAAatPos != aaToReplace){
    stop("AA to mutate different from original")
  }
  substr(aaSeqToMutate, aaReplacePos, aaReplacePos) <- aaReplacement
  results[i,]$mutant_aliphaticIndex <- aIndex(aaSeqToMutate)
  results[i,]$mutant_bomanIndex <- boman(aaSeqToMutate)
  results[i,]$mutant_charge <- charge(aaSeqToMutate)
  results[i,]$mutant_hydrophobicMoment <- hmoment(aaSeqToMutate)
  results[i,]$mutant_hydrophobicity <- hydrophobicity(aaSeqToMutate)
  results[i,]$mutant_instabilityIndex <- instaIndex(aaSeqToMutate)
  results[i,]$mutant_molWeight <- mw(aaSeqToMutate)
  results[i,]$mutant_massOverCharge <- mz(aaSeqToMutate)
  results[i,]$mutant_isoElecPoint <- pI(aaSeqToMutate)
  results[i,]$delta_aliphaticIndex <- results[i,]$mutant_aliphaticIndex - genePeptideProp$WT_aliphaticIndex
  results[i,]$delta_bomanIndex <- results[i,]$mutant_bomanIndex - genePeptideProp$WT_bomanIndex
  results[i,]$delta_charge <- results[i,]$mutant_charge - genePeptideProp$WT_charge
  results[i,]$delta_hydrophobicMoment <- results[i,]$mutant_hydrophobicMoment - genePeptideProp$WT_hydrophobicMoment
  results[i,]$delta_hydrophobicity <-  results[i,]$mutant_hydrophobicity - genePeptideProp$WT_hydrophobicity
  results[i,]$delta_instabilityIndex <- results[i,]$mutant_instabilityIndex - genePeptideProp$WT_instabilityIndex
  results[i,]$delta_molWeight <- results[i,]$mutant_molWeight - genePeptideProp$WT_molWeight
  results[i,]$delta_massOverCharge <- results[i,]$mutant_massOverCharge - genePeptideProp$WT_massOverCharge
  results[i,]$delta_isoElecPoint <- results[i,]$mutant_isoElecPoint - genePeptideProp$WT_isoElecPoint
  results[i,]$mutant_aaSeq <- aaSeqToMutate
}


###############################################################################
# Rename and move columns, add WT properties, and calculate mutant properties #
###############################################################################
# For clarity, rename mutant folding results to 'delta_x'
results <- results %>% rename(delta_total.energy = total.energy,
                              delta_Backbone.Hbond = Backbone.Hbond,
                              delta_Sidechain.Hbond = Sidechain.Hbond,
                              delta_Van.der.Waals = Van.der.Waals,
                              delta_Electrostatics = Electrostatics,
                              delta_Solvation.Polar = Solvation.Polar,
                              delta_Solvation.Hydrophobic = Solvation.Hydrophobic,
                              delta_Van.der.Waals.clashes = Van.der.Waals.clashes,
                              delta_entropy.sidechain = entropy.sidechain,
                              delta_entropy.mainchain = entropy.mainchain,
                              delta_sloop_entropy = sloop_entropy,
                              delta_mloop_entropy = mloop_entropy,
                              delta_cis_bond = cis_bond,
                              delta_torsional.clash = torsional.clash,
                              delta_backbone.clash = backbone.clash,
                              delta_helix.dipole = helix.dipole,
                              delta_water.bridge = water.bridge,
                              delta_disulfide = disulfide,
                              delta_electrostatic.kon = electrostatic.kon,
                              delta_partial.covalent.bonds = partial.covalent.bonds,
                              delta_energy.Ionisation = energy.Ionisation,
                              delta_Entropy.Complex = Entropy.Complex
                              )
# Add WT gene/protein properties to the variants
results <- merge(results, peptidePropPerGene, by.x = "gene", by.y = "WT_gene")
results$mutant_total.energy <- results$WT_total.energy + results$delta_total.energy
results$mutant_Backbone.Hbond <- results$WT_Backbone.Hbond + results$delta_Backbone.Hbond
results$mutant_Sidechain.Hbond <- results$WT_Sidechain.Hbond + results$delta_Sidechain.Hbond
results$mutant_Van.der.Waals <- results$WT_Van.der.Waals + results$delta_Van.der.Waals
results$mutant_Electrostatics <- results$WT_Electrostatics + results$delta_Electrostatics
results$mutant_Solvation.Polar <- results$WT_Solvation.Polar + results$delta_Solvation.Polar
results$mutant_Solvation.Hydrophobic <- results$WT_Solvation.Hydrophobic + results$delta_Solvation.Hydrophobic
results$mutant_Van.der.Waals.clashes <- results$WT_Van.der.Waals.clashes + results$delta_Van.der.Waals.clashes
results$mutant_entropy.sidechain <- results$WT_entropy.sidechain + results$delta_entropy.sidechain
results$mutant_entropy.mainchain <- results$WT_entropy.mainchain + results$delta_entropy.mainchain
results$mutant_sloop_entropy <- results$WT_sloop_entropy + results$delta_sloop_entropy
results$mutant_mloop_entropy <- results$WT_mloop_entropy + results$delta_mloop_entropy
results$mutant_cis_bond <- results$WT_cis_bond + results$delta_cis_bond
results$mutant_torsional.clash <- results$WT_torsional.clash + results$delta_torsional.clash
results$mutant_backbone.clash <- results$WT_backbone.clash + results$delta_backbone.clash
results$mutant_helix.dipole <- results$WT_helix.dipole + results$delta_helix.dipole
results$mutant_water.bridge <- results$WT_water.bridge + results$delta_water.bridge
results$mutant_disulfide <- results$WT_disulfide + results$delta_disulfide
results$mutant_electrostatic.kon <- results$WT_electrostatic.kon + results$delta_electrostatic.kon
results$mutant_partial.covalent.bonds <- results$WT_partial.covalent.bonds + results$delta_partial.covalent.bonds
results$mutant_energy.Ionisation <- results$WT_energy.Ionisation + results$delta_energy.Ionisation
results$mutant_Entropy.Complex <- results$WT_Entropy.Complex + results$delta_Entropy.Complex
# Remove Pdb file names, they are not useful
results$Pdb <- NULL
results$WT_Pdb <- NULL 
# Rename columns for clarity
results <- results %>% rename(dna_variant_chrom = chrom,
                              dna_variant_pos = pos,
                              dna_variant_ref = ref,
                              dna_variant_alt = alt,
                              dna_variant_assembly = assembly,
                              ann_classificationVKGL = classificationVKGL,
                              UniProtID = WT_UniProtKB.Swiss.Prot.ID,
                              TranscriptID = WT_Transcript.stable.ID,
                              ann_proteinLocalization = WT_proteinLocalization,
                              ann_proteinIschaperoned = WT_chaperoned,
                              ann_WT_energyPerKDa = WT_energyPerKDa,
                              ann_mutant_energy_SD = SD,
                              delta_aaSeq = protChange)
# Sort columns alphabetically and move selected ones to the front
results <- results[,order(colnames(results))]
results <- results %>% relocate(gene, TranscriptID, UniProtID, dna_variant_chrom, dna_variant_pos, dna_variant_ref, dna_variant_alt, dna_variant_assembly, ann_classificationVKGL)
colnames(results)


#####################
# Add QSAR features #
#####################
# Better rownames to results, makes it easier to merge later
#rownames(results) <- paste0(results$gene, "/", results$UniProtID, ":", results$delta_aaSeq)
mutant_qsarFt <- data.frame()
for(i in 1:nrow(results))
{
  mutantAAseq <- results[i,'mutant_aaSeq']
  qsarRes <- peptoolkit::extract_features_QSAR(nchar(mutantAAseq), custom.list = TRUE, PeList = c(mutantAAseq))
  mutant_qsarFt <- rbind(mutant_qsarFt, qsarRes)
}
# repeat for WT
WT_qsarFt <- data.frame()
for(i in 1:nrow(results))
{
  WTAAseq <- results[i,'WT_aaSeq']
  qsarRes <- peptoolkit::extract_features_QSAR(nchar(WTAAseq), custom.list = TRUE, PeList = c(WTAAseq))
  WT_qsarFt <- rbind(WT_qsarFt, qsarRes)
}
# Write/read to buffer intermediate data
#write.csv(mutant_qsarFt, paste(rootDir, "data", "mutant_qsarFt.csv", sep="/"), row.names = FALSE, quote = FALSE)
#write.csv(WT_qsarFt, paste(rootDir, "data", "WT_qsarFt.csv", sep="/"), row.names = FALSE, quote = FALSE)
#mutant_qsarFt <- read.csv(paste(rootDir, "data", "mutant_qsarFt.csv", sep="/"))
#WT_qsarFt <- read.csv(paste(rootDir, "data", "WT_qsarFt.csv", sep="/"))
mutant_qsarFt$Sequence <- NULL # remove sequence, identical to mutant_aaSeq
WT_qsarFt$Sequence <- NULL # remove sequence, identical to WT_aaSeq
delta_qsarFt <- mutant_qsarFt-WT_qsarFt # compute deltas
# Rename columns
colnames(mutant_qsarFt) <- paste("mutant", colnames(mutant_qsarFt), sep = "_")
colnames(WT_qsarFt) <- paste("WT", colnames(WT_qsarFt), sep = "_")
colnames(delta_qsarFt) <- paste("delta", colnames(delta_qsarFt), sep = "_")
# Combine with results
results <- cbind(results, delta_qsarFt, mutant_qsarFt, WT_qsarFt)


#######################################################
# Pseudo-Amino Acid Composition (PseAAC) and.         #
# Amphiphilic Pseudo Amino Acid Composition (APseAAC) #
#######################################################
res_APseAAC <- data.frame()
for(i in 1:nrow(results))
{
  cat(paste0("working on row ", i, "\n"))
  wtAAseq <- results[i,'WT_aaSeq']
  mutantAAseq <- results[i,'mutant_aaSeq']
  wtPAAC <- as.data.frame(t(protr::extractPAAC(wtAAseq)))
  wtAPAAC <- as.data.frame(t(protr::extractAPAAC(wtAAseq)))
  mutantPAAC <- as.data.frame(t(protr::extractPAAC(mutantAAseq)))
  mutantAPAAC <- as.data.frame(t(protr::extractAPAAC(mutantAAseq)))
  deltaPAAC <- mutantPAAC-wtPAAC
  deltaAPAAC <- mutantAPAAC-wtAPAAC
  colnames(wtPAAC) <- paste("WT", colnames(wtPAAC), sep = "_")
  colnames(wtAPAAC) <- paste("WT", colnames(wtAPAAC), sep = "_")
  colnames(mutantPAAC) <- paste("mutant", colnames(mutantPAAC), sep = "_")
  colnames(mutantAPAAC) <- paste("mutant", colnames(mutantAPAAC), sep = "_")
  colnames(deltaPAAC) <- paste("delta", colnames(deltaPAAC), sep = "_")
  colnames(deltaAPAAC) <- paste("delta", colnames(deltaAPAAC), sep = "_")
  allCols <- cbind(wtPAAC, wtAPAAC, mutantPAAC, mutantAPAAC, deltaPAAC, deltaAPAAC)
  res_APseAAC <- rbind(res_APseAAC, allCols)
}
write.csv(res_APseAAC, paste(rootDir, "data", "res_APseAAC.csv", sep="/"), row.names = FALSE, quote = FALSE)
results <- cbind(results, res_APseAAC)


#######################################################
# Write/read complete collation of variant level data #
#######################################################
freeze3 <- paste(rootDir, "data", "freeze3.csv.gz", sep="/")
#write.csv.gz(results, freeze3, row.names = FALSE, quote = FALSE)
results <- read.csv(freeze3)
# NOTE: freeze3 data was superseded by freeze4 and therefore left out


######################################
# Add AlphaMissense am_pathogenicity #
######################################
# see: CombineWithAlhaMissense.R
# This results in the creation of freeze4.csv.gz


######################################
# Functional annotations
######################################
# After this:
# --> main-ligand-binding-sites.R
#     Runs P2Rank and adds <mutations>_<protein>_v4_Repair.pdb_predictions.csv
#     Contains predicted ligand binding pockets
# --> main-dna-rna-binding.R
#     Runs GLM-Score and adds DNA_RNA_interaction_terms.csv
#     Contains simulated binding energy to actual DNA and RNA molecules
# --> main-geo-net.R
#     Runs GeoNet and adds CombinedGeoNetPredictions.csv
#     Contains predicted protein, DNA and RNA binding sites
#
# If all is complete, based on freeze4 dataset:
# 2596 wildtypes + 23417 mutations = 26013 files of each:
#     P2Rank:     find . -name "*_v4_Repair.pdb_predictions.csv" | wc -l
#     GLM-Score:  find . -name "DNA_RNA_interaction_terms.csv" | wc -l
#     GeoNet:     find . -name "CombinedGeoNetPredictions.csv" | wc -l
#
# make-freeze5.R combines these functional features and remove non-functional ones
# Resulting in freeze5.csv.gz
#
