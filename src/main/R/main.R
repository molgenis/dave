######################
# Load used packages #
######################
library(R.utils)   # for 'gunzip', 'mkdirs'
library(Rpdb)      # to load PDB files
library(dplyr)     # to remove duplicate rows
library(Peptides)  # protein annotations


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
secr <- read.table(file=paste(rootDir, "data", "protein-atlas-secreted-genenames-mane-uniprot-withvariants.txt", sep="/"), sep = '\t',header = TRUE)
secr$protType <- "secreted"
intr <- read.table(file=paste(rootDir, "data", "protein-atlas-intracellular-genenames-mane-uniprot-random2000-withvariants.txt", sep="/"), sep = '\t',header = TRUE)
intr$protType <- "intracellular"
memb <- read.table(file=paste(rootDir, "data", "protein-atlas-membrane-genenames-mane-uniprot-random2000-withvariants.txt", sep="/"), sep = '\t',header = TRUE)
memb$protType <- "membrane"
allg <- rbind(secr, intr, memb)
# Selected gene list to work on
selectedGenes <- allg


##########################################
# Derived paths of directories and files #
##########################################
dataGenesDir <- paste(rootDir, "data", "genes", sep="/")
source(paste(rootDir, "src", "main", "R", "aa3to1.R", sep="/"))


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
      cat("  already tried before but failed, breaking...\n")
      break
    }
    file.copy(from = repPDB_fullLoc, to = mutationDir)
    write(paste(vkgl[j, "ProtChange"], ";", sep=""), file = "individual_list.txt")
    state <- system(paste(foldxExec, " --command=BuildModel --mutant-file=individual_list.txt --numberOfRuns=5 --pdb=", repPDB_fileName, sep=""), intern = TRUE)
    if(any(grepl("Specified residue not found", state)))
      {
        write(state, file = "exception.txt")
        break
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
    mutationToCheck <- vkgl[1,]$ProtChange
    mutationResultsDir <- paste(specificGeneDir, "folding-results", mutationToCheck, sep="/")
    avgEnergy <- list.files(mutationResultsDir, pattern="Average")
    if(length(avgEnergy) > 0){
      cat(paste("Gene", geneName, "has at least 1 successful mutation:", mutationToCheck,"so doing quick refold to get WT DG...\n", sep=" "))
      tmpDir <- paste(specificGeneDir, "tmp", sep="/")
      mkdirs(tmpDir)
      setwd(tmpDir)
      file.copy(from = paste(specificGeneDir, pdbFile, sep="/"), to = tmpDir)
      write(paste(mutationToCheck, ";", sep=""), file = "individual_list.txt")
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
        stop(paste0("Something went wrong with gene ", geneName, ", mutation ", + mutationToCheck))
      }
    }else{
      cat(paste("No average energy file, mutation failed ? Skipping", geneName, "\n", sep=" "))
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
    geneInfo$wtDG <- wtInfo[1,]$total.energy
  }else{
    geneInfo$wtDG <- NA
  }
  variants <- read.table(file=paste(specificGeneDir, vkglProtVarFileName, sep="/"), sep = '\t', header = TRUE)
  foldingResultsDir <- paste(specificGeneDir, "folding-results", sep="/")
  for(j in seq_along(variants))
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
    #result$vkglGeneNameAnnotation <- variants[j, "Gene"]
    result$gene <- geneName
    result$protChange <- variants[j, "ProtChange"]
    result$classificationVKGL <- variants[j, "Classification"]
    result$transcript <- geneInfo$Transcript.stable.ID
    result$uniprot <- geneInfo$UniProtKB.Swiss.Prot.ID
    result$protType <- geneInfo$protType
    result$wtDG <- geneInfo$wtDG
    results <- rbind(results, result)
  }
}


################################################
# Retrieve genes that interact with chaperones #
################################################
chapInteractingGenesLoc <- paste(rootDir, "data", "chaperones", "interacting-with-chaperones-genenames-merged-with-uniprotmapped.txt", sep="/")
chapInteractingGenes <- read.table(file=chapInteractingGenesLoc, sep = '\t',header = TRUE)
results$chaperoned <- as.factor(results$gene %in% chapInteractingGenes$Gene.name)
#alternative way to add annotation? something like
#result$aaa <- as.factor(selectedGenes[selectedGenes$Gene.name==results$gene, "protType"])

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
    
    cat(paste("calculations done\n", sep=" "))

  }else{
    cat(paste("No PDB file for gene", geneName, "\n", sep=" "))
  }
}
# Persist these results for quick loading later
peptidePropPerGeneFile <- paste(rootDir, "data", "peptidePropPerGene.csv", sep="/")
#write.csv(peptidePropPerGene, peptidePropPerGeneFile, row.names = FALSE, quote = FALSE)
peptidePropPerGene <- read.csv(peptidePropPerGeneFile)


######3
# For all succesfully folded mutations, introduce AA change in AA sequence and calculate additional properties
#######
for(i in seq_along(results))
{
 # todo
}

############################################################################################################
# Write/read data freeze1 based on all secreted and a random selection of membrane and intracellular genes #
############################################################################################################
freeze1 <- paste(rootDir, "data", "freeze1.csv", sep="/")
write.csv(results, freeze1, row.names = FALSE, quote = FALSE)
results <- read.csv(freeze1)
