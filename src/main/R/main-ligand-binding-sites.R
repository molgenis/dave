######################
# Load used packages #
######################
library(R.utils)   # for 'gunzip', 'mkdirs'

#######################
# Adjustable settings #
#######################
# Example for Unix
rootDir <- "/Users/joeri/git/dave" # root directory that contains README.md, data/, img/, out/, src/, etc.
vkglProtVarFileName <- "VKGL_apr2024_protForFolding.tsv" # loading processed VKGL protein variants from: /data/genes/{gene}/{vkglProtVarFile}
foldxExec <- "/Applications/FoldX/5/foldx5MacStd/foldx_20241231" # exact path to the FoldX executable
p2rankExec <- "/Applications/p2rank_2.5/prank"
# Example for Windows
rootDir <- "D:/github/dave"
vkglProtVarFileName <- "VKGL_apr2024_protForFolding.tsv"
foldxExec <- "C:/\"Program Files\"/FoldX/foldx_20241231.exe"
java <- r"(C:\Program Files\OpenLogic\jdk-22.0.2.9-hotspot\bin\java.exe)"
prank <- r"(D:/p2rank_2.5/bin/p2rank.jar;D:/p2rank_2.5/bin/lib/*)"

# Assuming data was produced by main.R and enriched by CombineWithAlhaMissense.R
# Load the freeze4 data and assign meaningful row names
freeze <- paste(rootDir, "data", "freeze4.csv.gz", sep="/")
results <- read.csv(freeze)
rownames(results) <- paste0(results$gene, "/", results$UniProtID, ":", results$delta_aaSeq)
dataGenesDir <- paste(rootDir, "data", "genes", sep="/")
source(paste(rootDir, "src", "main", "R", "helper-functions.R", sep="/"))

# Iterate over genes and then over mutations
setwd(dataGenesDir)
succesfulGenes <- unique(results$gene)

# Possible gene lists to work on
secr <- read.table(file=paste(rootDir, "data", "protein-atlas-secreted-geneIDs-mane-uniprot-withvariants.txt", sep="/"), sep = '\t',header = TRUE)
intr <- read.table(file=paste(rootDir, "data", "protein-atlas-intracellular-geneIDs-mane-uniprot-random2000-withvariants.txt", sep="/"), sep = '\t',header = TRUE)
memb <- read.table(file=paste(rootDir, "data", "protein-atlas-membrane-geneIDs-mane-uniprot-random2000-withvariants.txt", sep="/"), sep = '\t',header = TRUE)
succesfulGenesSub <- Reduce(intersect, list(secr$Gene.name, succesfulGenes))

# Iterate over selection of genes
for(i in seq_along(succesfulGenesSub))
{
  #i <- 1 # DEBUG/DEV
  geneName <- succesfulGenesSub[i]
  cat(paste("Working on gene:", geneName, "\n", sep=" "))
  specificGeneDir <- paste(dataGenesDir, geneName, sep="/")
  setwd(specificGeneDir)
  pdbFile <- list.files(pattern="*_Repair.pdb$")
  if(length(pdbFile) == 0){
    stop(paste("No PDB file for gene", geneName, "\n", sep=" "))
  }
  if(length(list.files(specificGeneDir, pattern="Repair.pdb_predictions.csv")) == 0){
    # ligand binding sites for WT structure, keep only _predictions file
    tmpDir <- paste(specificGeneDir, "tmp", sep="/")
    mkdirs(tmpDir)
    setwd(tmpDir)
    file.copy(from = paste(specificGeneDir, pdbFile, sep="/"), to = tmpDir)
    prArgs <- paste0(" predict -c alphafold -f ", pdbFile ," -o ",tmpDir)
    
    if(.Platform$OS.type == "windows")
    {
      system2(command = java, args = c("-Xmx2048m", paste0("-cp ", prank), "cz.siret.prank.program.Main", prArgs) )
    }else{
      system(paste(p2rankExec, prArgs ,sep=""), intern = TRUE)
    }
    
    wtPred <- paste0(pdbFile,"_predictions.csv")
    file.copy(from = paste(tmpDir, wtPred, sep="/"), to = specificGeneDir)
    if(!file.exists(paste(specificGeneDir,wtPred,sep="/"))){
      stop(paste("Missing wild-type p2rank output:", wtPred, "\n", sep=" "))
    }
      setwd(specificGeneDir)
      unlink(tmpDir, recursive = TRUE)
  }else{
    cat("Wild-type *Repair.pdb_predictions.csv file already present, skipping..\n")
  }
  
  # Iterate over variants, mutate PDB and predict ligand binding sites
  variants <- read.table(file=paste(specificGeneDir, vkglProtVarFileName, sep="/"), sep = '\t', header = TRUE, colClasses = c("character", "character", "numeric", "character", "character", "character", "character", "numeric", "character"))
  foldingResultsDir <- paste(specificGeneDir, "folding-results", sep="/")
  for(j in 1:nrow(variants))
  {
    #j <- 2 # DEBUG/DEV
    mutation <- variants$ProtChange[j]
    cat(paste("Working on ", mutation, " (gene ",geneName,", mutation ", j, " of ", nrow(variants), ")\n", sep=""))
    mutationDir <- paste(foldingResultsDir, mutation, sep="/")
    if(!dir.exists(mutationDir))
    {
      stop(paste("No mutation dir", mutationDir, "\n", sep=" "))
    }
    if(!length(list.files(mutationDir, pattern="Repair.pdb_predictions.csv")) == 0){
      cat("Mutant *Repair.pdb_predictions.csv file already present, skipping..\n")
      next
    }
    if(!length(list.files(mutationDir, pattern="exception.txt")) == 0){
      cat("Exception present, skipping..\n")
      next
    }
    tmpDir <- paste(mutationDir, "tmp", sep="/")
    mkdirs(tmpDir)
    setwd(tmpDir)
    file.copy(from = paste(specificGeneDir, pdbFile, sep="/"), to = tmpDir)
    system(paste(foldxExec, " --command=PositionScan --pdb=", pdbFile, " --out-pdb=true --positions=",mutation,sep=""), intern = TRUE)
    mutantAA1 <- substr(mutation, nchar(mutation), nchar(mutation))
    position <- substr(mutation, 3, nchar(mutation)-1)
    mutantAA3 <- aa1to3(mutantAA1)
    mutantPDB <- paste0(mutantAA3, position, "_", pdbFile)
    
    prArgs <- paste0(" predict -c alphafold -f ", mutantPDB ," -o ", tmpDir)
    
    if(.Platform$OS.type == "windows")
    {
      system2(command = java, args = c("-Xmx2048m", paste0("-cp ", prank), "cz.siret.prank.program.Main", prArgs) )
    }else{
      system(paste(p2rankExec, prArgs, sep=""), intern = TRUE)
    }

    mutantLBSPred <- paste0(mutantAA3, position, "_", pdbFile, "_predictions.csv")
    file.copy(from = paste(tmpDir, mutantLBSPred, sep="/"), to = mutationDir)
    if(!file.exists(paste(mutationDir, mutantLBSPred, sep="/"))){
      stop(paste("Missing mutant p2rank output:", mutantLBSPred, "\n", sep=" "))
    }
    setwd(mutationDir)
    unlink(tmpDir, recursive = TRUE)
  }
}


###########################################################
# Run a particular mutation manually, keeping all outputs #
###########################################################
# --> This is a "hardcoded" example on how to run this on cmdline!
#
# cd /Users/joeri/git/dave/data/genes/SLCO2A1/folding-results/GA554R
# mkdir tmp
# cp ../../AF-Q92959-F1-model_v4_Repair.pdb tmp 
# cd tmp
# /Applications/FoldX/5/foldx_20251231 --command=PositionScan --pdb=AF-Q92959-F1-model_v4_Repair.pdb --out-pdb=true --positions=GA554R
# /Applications/p2rank_2.5/prank predict -c alphafold -f ARG554_AF-Q92959-F1-model_v4_Repair.pdb -o .
# 
# 
# if error "Unsupported class file major version 68", lower the Java version:
# /usr/libexec/java_home -V
# export JAVA_HOME=`/usr/libexec/java_home -v 17`
#
