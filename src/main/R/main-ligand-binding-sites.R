######################
# Load used packages #
######################
library(R.utils)   # for 'gunzip', 'mkdirs'

#######################
# Adjustable settings #
#######################
# Example for Unix
rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability" # root directory that contains README.md, data/, img/, out/, src/, etc.
vkglProtVarFileName <- "VKGL_apr2024_protForFolding.tsv" # loading processed VKGL protein variants from: /data/genes/{gene}/{vkglProtVarFile}
foldxExec <- "/Applications/FoldX/5/foldx5MacStd/foldx_20241231" # exact path to the FoldX executable
p2rankExec <- "/Applications/p2rank_2.5/prank"
# Example for Windows
rootDir <- "D:/github/vkgl-secretome-protein-stability"
vkglProtVarFileName <- "VKGL_apr2024_protForFolding.tsv"
foldxExec <- "C:/\"Program Files\"/FoldX/foldx_20241231.exe"
p2rankExec <- "C:/\"Program Files\"/p2rank_2.5/prank"

# Assuming data was produced by main.R and enriched by CombineWithAlhaMissense.R
# Load the freeze4 data and assign meaningful row names
freeze <- paste(rootDir, "data", "freeze4.csv.gz", sep="/")
results <- read.csv(freeze)
rownames(results) <- paste0(results$gene, "/", results$UniProtID, ":", results$delta_aaSeq)
dataGenesDir <- paste(rootDir, "data", "genes", sep="/")
source(paste(rootDir, "src", "main", "R", "aa3to1.R", sep="/"))

# Iterate over genes and then over mutations
setwd(dataGenesDir)
succesfulGenes <- unique(results$gene)

for(i in seq_along(succesfulGenes))
{
  i <- 1 # DEBUG/DEV
  geneName <- succesfulGenes[i]
  
  cat(paste("working on gene:", geneName, "\n", sep=" "))
  specificGeneDir <- paste(dataGenesDir, geneName, sep="/")
  setwd(specificGeneDir)
  pdbFile <- list.files(pattern="*_Repair.pdb")
  if(length(pdbFile) == 0){
    stop(paste("No PDB file for gene", geneName, "\n", sep=" "))
  }
  
  # predict on WT structure
  tmpDir <- paste(specificGeneDir, "tmp", sep="/")
  mkdirs(tmpDir)
  setwd(tmpDir)
  file.copy(from = paste(specificGeneDir, pdbFile, sep="/"), to = tmpDir)
  system(paste(p2rankExec, " predict -c alphafold -f ", pdbFile ," -o ",tmpDir,sep=""), intern = TRUE)
  
  # TAR GZ results
  
  
  
  variants <- read.table(file=paste(specificGeneDir, vkglProtVarFileName, sep="/"), sep = '\t', header = TRUE, colClasses = c("character", "character", "numeric", "character", "character", "character", "character", "numeric", "character"))
  foldingResultsDir <- paste(specificGeneDir, "folding-results", sep="/")
  
  for(j in 1:nrow(variants))
  {
    j <- 2 # DEBUG/DEV
    mutation <- variants$ProtChange[j]
    cat(paste("working on mutant:", geneName, mutation, "\n", sep=" "))
    
    mutationDir <- paste(foldingResultsDir, mutation, sep="/")
    if(!dir.exists(mutationDir))
    {
      stop(paste("No mutation dir", mutationDir, "\n", sep=" "))
    }
    if(length(list.files(mutationDir, pattern="*.fxout")) == 0){
      next
    }
    
    tmpDir <- paste(mutationDir, "tmp", sep="/")
    mkdirs(tmpDir)
    setwd(tmpDir)
    file.copy(from = paste(specificGeneDir, pdbFile, sep="/"), to = tmpDir)
    
    system(paste(foldxExec, " --command=PositionScan --pdb=",pdbFile,"  --out-pdb=true --positions=",mutation,sep=""), intern = TRUE)
    
    wtAA <- substr(mutation, 1, 1)
    muAA <- substr(mutation, nchar(mutation), nchar(mutation))
    posi <- substr(mutation, 3, nchar(mutation)-1)
    wtAA3 <- aa1to3(wtAA)
    muAA3 <- aa1to3(muAA)
    
    wtPDB <- paste0(wtAA3,posi,"_",pdbFile)
    muPDB <- paste0(muAA3,posi,"_",pdbFile)

    system(paste(p2rankExec, " predict -c alphafold -f ", muPDB ," -o ",tmpDir,sep=""), intern = TRUE)
    
    wtPred <- paste0(wtAA3,posi,"_",pdbFile,"_predictions.csv")
    wtResi <- paste0(wtAA3,posi,"_",pdbFile,"_residues.csv")
    muPred <- paste0(muAA3,posi,"_",pdbFile,"_predictions.csv")
    muResi <- paste0(muAA3,posi,"_",pdbFile,"_residues.csv")
    
    if(!file.exists(wtPred) | !file.exists(wtResi) | !file.exists(muPred) | !file.exists(muResi)){
      stop(paste("Missing p2rank output for ", mutationDir, "\n", sep=" "))
    }
    
    filesToArchive <- c(wtPred, wtResi, muPred, muResi)
    tarFile <- paste0(mutation,".tar")
    tar(tarfile = tarFile, files = filesToArchive)
    gzip(tarFile, overwrite = TRUE)
    tarGzFile <- paste0(tarFile,".gz")
    file.copy(from = paste(tmpDir, tarGzFile, sep="/"), to = mutationDir)
    unlink(tmpDir, recursive = TRUE)
  }
}

