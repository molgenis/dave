######################
# Load used packages #
######################
library(R.utils)    # for 'gunzip', 'mkdirs'

#######################
# Adjustable settings #
#######################
# Example for Unix
rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability" # root directory that contains README.md, data/, img/, out/, src/, etc.
vkglProtVarFileName <- "VKGL_apr2024_protForFolding.tsv" # loading processed VKGL protein variants from: /data/genes/{gene}/{vkglProtVarFile}
foldxExec <- "/Applications/FoldX/5/foldx5MacStd/foldx_20241231" # exact path to the FoldX executable
glmScoreDir <- "/Users/joeri/git/GLM-Score"
glmScoreExe <- "GLM-Score"
# Example for Windows
rootDir <- "D:/github/vkgl-secretome-protein-stability"
vkglProtVarFileName <- "VKGL_apr2024_protForFolding.tsv"
foldxExec <- "C:/\"Program Files\"/FoldX/foldx_20241231.exe"
glmScoreDir <- "D:/github/GLM-Score"
glmScoreExe <- "GLM-Score"


# Assuming data was produced by main.R and enriched by CombineWithAlhaMissense.R
# Load the freeze4 data and assign meaningful row names
freeze <- paste(rootDir, "data", "freeze4.csv.gz", sep="/")
results <- read.csv(freeze)
rownames(results) <- paste0(results$gene, "/", results$UniProtID, ":", results$delta_aaSeq)
dataGenesDir <- paste(rootDir, "data", "genes", sep="/")
dataStructDir <- paste(rootDir, "data", "structures", sep="/")
glmScoreDirExe <- paste(glmScoreDir, glmScoreExe, sep="/")
resultFileName <- "DNA_RNA_interaction_terms.csv"
source(paste(rootDir, "src", "main", "R", "aa3to1.R", sep="/"))
source(paste(rootDir, "src", "main", "R", "combineDNAandRNAresults.R", sep="/"))


# Iterate over genes and then over mutations
setwd(dataGenesDir)
succesfulGenes <- unique(results$gene)

# Possible gene lists to work on
secr <- read.table(file=paste(rootDir, "data", "protein-atlas-secreted-geneIDs-mane-uniprot-withvariants.txt", sep="/"), sep = '\t',header = TRUE)
intr <- read.table(file=paste(rootDir, "data", "protein-atlas-intracellular-geneIDs-mane-uniprot-random2000-withvariants.txt", sep="/"), sep = '\t',header = TRUE)
memb <- read.table(file=paste(rootDir, "data", "protein-atlas-membrane-geneIDs-mane-uniprot-random2000-withvariants.txt", sep="/"), sep = '\t',header = TRUE)
succesfulGenesSub <- Reduce(intersect, list(secr$Gene.name, succesfulGenes))
results <- NULL # unload results to save memory

# DNA, in PDB and MOL2 format
file.copy(from = paste(dataStructDir, "1bna.pdb", sep="/"), to = glmScoreDir)
file.copy(from = paste(dataStructDir, "1bna.mol2", sep="/"), to = glmScoreDir)
# RNA, in PDB and MOL2 format
file.copy(from = paste(dataStructDir, "4jrd-A.pdb", sep="/"), to = glmScoreDir)
file.copy(from = paste(dataStructDir, "4jrd-A.mol2", sep="/"), to = glmScoreDir)

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
  if(!length(list.files(specificGeneDir, pattern=resultFileName)) == 0){
    cat("Wild-type combined result file already present, skipping..\n")
    
  }else{
    file.copy(from = paste(specificGeneDir, pdbFile, sep="/"), to = glmScoreDir)
    setwd(glmScoreDir)
    dna_pKd <- system(paste(glmScoreDirExe, paste0(pdbFile, " 1bna.pdb  1bna.mol2 DNA")), intern = TRUE)
    rna_pKd <- system(paste(glmScoreDirExe, paste0(pdbFile, " 4jrd-A.pdb 4jrd-A.mol2 RNA")), intern = TRUE)
    # Grab interaction terms for both and combine into one file. From https://github.com/Klab-Bioinfo-Tools/GLM-Score doc:
    # Total hydrophobic contact score (V2), Van der Waals interactions (V3), side chain rotation (V4), hydrogen bonding (V5), assessible to solvent area ASA of protein (V6) and ligand (V7), repulsive interactions (V18), london disperson forces (V19), contact hydrophobicity (V20), total hydrophobicity (V21), contact surface tension (V22), total surface tension (V23)
    transpose_and_label("1bna.interaction_terms.txt", "4jrd-A.interaction_terms.txt", dna_pKd, rna_pKd, resultFileName)
    file.copy(from = resultFileName, to = specificGeneDir)
    if(!file.exists(paste(specificGeneDir, resultFileName, sep="/"))){ stop("Copy of combined result file failed") }
    unlink("1bna.interaction_terms.txt")
    unlink("4jrd-A.interaction_terms.txt")
    unlink(resultFileName)
    if(file.exists("1bna.interaction_terms.txt")){ stop("Old 1bna.interaction_terms.txt still present") }
    if(file.exists("4jrd-A.interaction_terms.txt")){ stop("Old 4jrd-A.interaction_terms.txt still present") }
    if(file.exists(resultFileName)){ stop("Old combined result file still present") }
  }
  
  # Iterate over variants, mutate PDB and predict ligand binding sites
  variants <- read.table(file=paste(specificGeneDir, vkglProtVarFileName, sep="/"), sep = '\t', header = TRUE, colClasses = c("character", "character", "numeric", "character", "character", "character", "character", "numeric", "character"))
  foldingResultsDir <- paste(specificGeneDir, "folding-results", sep="/")
  for(j in 1:nrow(variants))
  {
    #j <- 1 # DEBUG/DEV
    mutation <- variants$ProtChange[j]
    cat(paste("Working on ", mutation, " (gene ",geneName,", mutation ", j, " of ", nrow(variants), ")\n", sep=""))
    mutationDir <- paste(foldingResultsDir, mutation, sep="/")
    if(!dir.exists(mutationDir))
    {
      stop(paste("No mutation dir", mutationDir, "\n", sep=" "))
    }
    if(!length(list.files(mutationDir, pattern=resultFileName)) == 0){
      cat("Mutant combined result file already present, skipping..\n")
      next
    }
    if(!length(list.files(mutationDir, pattern="exception.txt")) == 0){
      cat("Exception present, skipping..\n")
      next
    }
    
    tmpDir <- paste(mutationDir, "tmp-dna-rna", sep="/")
    mkdirs(tmpDir)
    setwd(tmpDir)
    file.copy(from = paste(specificGeneDir, pdbFile, sep="/"), to = tmpDir)
    
    system(paste(foldxExec, " --command=PositionScan --pdb=", pdbFile, " --out-pdb=true --positions=", mutation,sep=""), intern = TRUE)
    mutantAA1 <- substr(mutation, nchar(mutation), nchar(mutation))
    position <- substr(mutation, 3, nchar(mutation)-1)
    mutantAA3 <- aa1to3(mutantAA1)
    mutantPDB <- paste0(mutantAA3, position, "_", pdbFile)
    
    if(!file.exists(mutantPDB)){
      Sys.sleep(1)
      system(paste(foldxExec, " --command=PositionScan --pdb=", pdbFile, " --out-pdb=true --positions=", mutation,sep=""), intern = TRUE)
    }
    
    file.copy(from = paste(tmpDir, mutantPDB, sep="/"), to = glmScoreDir)
    setwd(glmScoreDir)
    dna_pKd <- system(paste(glmScoreDirExe, paste0(mutantPDB, " 1bna.pdb  1bna.mol2 DNA")), intern = TRUE)
    rna_pKd <- system(paste(glmScoreDirExe, paste0(mutantPDB, " 4jrd-A.pdb 4jrd-A.mol2 RNA")), intern = TRUE)
    transpose_and_label("1bna.interaction_terms.txt", "4jrd-A.interaction_terms.txt", dna_pKd, rna_pKd, resultFileName)
    file.copy(from = resultFileName, to = mutationDir)
    if(!file.exists(paste(mutationDir, resultFileName, sep="/"))){ stop("Copy of combined result file failed") }
    unlink("1bna.interaction_terms.txt")
    unlink("4jrd-A.interaction_terms.txt")
    unlink(resultFileName)
    if(file.exists("1bna.interaction_terms.txt")){ stop("Old 1bna.interaction_terms.txt still present") }
    if(file.exists("4jrd-A.interaction_terms.txt")){ stop("Old 4jrd-A.interaction_terms.txt still present") }
    if(file.exists(resultFileName)){ stop("Old combined result file still present") }
    unlink(tmpDir, recursive = TRUE)
  }
}


