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
pythonExe <- "/opt/homebrew/anaconda3/bin/python"
pdbFixerPath <- "/opt/homebrew/Caskroom/miniconda/base/bin/pdbfixer" # installed via: conda install -c conda-forge pdbfixer
geoNetScriptsDir <- "/Users/joeri/git/GeoNet/scripts"
geoNetExe <- "prediction.py"
nrOfCores <- 6
# Example for Windows
rootDir <- "D:/github/vkgl-secretome-protein-stability"
vkglProtVarFileName <- "VKGL_apr2024_protForFolding.tsv"
foldxExec <- "C:/\"Program Files\"/FoldX/foldx_20241231.exe"
geoNetScriptsDir <- "D:/github/GeoNet/scripts"
geoNetExe <- "prediction.py"

# Assuming data was produced by main.R and enriched by CombineWithAlhaMissense.R
# Load the freeze4 data and assign meaningful row names
freeze <- paste(rootDir, "data", "freeze4.csv.gz", sep="/")
results <- read.csv(freeze)
rownames(results) <- paste0(results$gene, "/", results$UniProtID, ":", results$delta_aaSeq)
dataGenesDir <- paste(rootDir, "data", "genes", sep="/")
source(paste(rootDir, "src", "main", "R", "aa3to1.R", sep="/"))
source(paste(rootDir, "src", "main", "R", "replace-line-in-file.R", sep="/"))
source(paste(rootDir, "src", "main", "R", "combine-last-columns.R", sep="/"))

# Iterate over genes and then over mutations
setwd(dataGenesDir)
succesfulGenes <- unique(results$gene)

# Possible gene lists to work on
secr <- read.table(file=paste(rootDir, "data", "protein-atlas-secreted-geneIDs-mane-uniprot-withvariants.txt", sep="/"), sep = '\t',header = TRUE)
intr <- read.table(file=paste(rootDir, "data", "protein-atlas-intracellular-geneIDs-mane-uniprot-random2000-withvariants.txt", sep="/"), sep = '\t',header = TRUE)
memb <- read.table(file=paste(rootDir, "data", "protein-atlas-membrane-geneIDs-mane-uniprot-random2000-withvariants.txt", sep="/"), sep = '\t',header = TRUE)
succesfulGenesSub <- Reduce(intersect, list(secr$Gene.name, succesfulGenes))
results <- NULL # unload results to save memory


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
  if(!length(list.files(specificGeneDir, pattern="CombinedGeoNetPredictions.csv")) == 0){
    cat("Wild-type CombinedGeoNetPredictions.csv file already present, skipping..\n")
  }else{
    tmpDir <- paste(specificGeneDir, "tmp-geonet", sep="/")
    mkdirs(tmpDir)
    setwd(tmpDir)
    unfixedPdb <-  paste0(tmpDir, "/", pdbFile, ".unfixed")
    file.copy(from = paste(specificGeneDir, pdbFile, sep="/"), to = unfixedPdb)
    tmpPDBPath <- paste(tmpDir, pdbFile, sep="/")
    system(paste(pdbFixerPath, " ", unfixedPdb, " --output ", pdbFile, sep=""), intern = TRUE)
    setwd(geoNetScriptsDir)
    geoNetExe_WT <- paste0(geneName,"_WT", geoNetExe)
    file.copy(from = geoNetExe, to = geoNetExe_WT)
    oldArgs <- "    args = parse_args()"
    newArgs_int <- paste0("    args = parse_args([\"--querypath\", \"",tmpDir,"\", \"--filename\", \"",pdbFile,"\", \"--chainid\", \"A\", \"--ligand\", \"P,DNA,RNA\", \"--cpu\", \"",nrOfCores,"\"])")
    ok <- replace_line_in_file(geoNetExe_WT, oldArgs, newArgs_int)
    if(!ok){ stop("Args replace failed") }
    system(paste(pythonExe, geoNetExe_WT, sep=" "), intern = TRUE)
    unlink(geoNetExe_WT)
    setwd(tmpDir)
    dnaFile <- paste0(tmpDir, "/", list.files(pattern = "/*DNA-binding_result.csv"))
    rnaFile <- paste0(tmpDir, "/", list.files(pattern = "/*RNA-binding_result.csv"))
    pFile <- paste0(tmpDir, "/", list.files(pattern = "/*P-binding_result.csv"))
    combinedFile <- "CombinedGeoNetPredictions.csv"
    combine_geonet_results(dnaFile, rnaFile, pFile, c("dnaProb", "dnaBin", "rnaProb", "rnaBin", "pProb", "pBin"), combinedFile)
    if(!file.exists(combinedFile)){ stop("Combined file for WT failed") }
    file.copy(from = combinedFile, to = specificGeneDir)
    if(!file.exists(paste(specificGeneDir, combinedFile, sep="/"))){ stop("Combined file for WT not in final location") }
    setwd(specificGeneDir)
    unlink(tmpDir, recursive = TRUE)
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
    if(!length(list.files(mutationDir, pattern="CombinedGeoNetPredictions")) == 0){
      cat("Mutant CombinedGeoNetPredictions file already present, skipping..\n")
      next
    }
    if(!length(list.files(mutationDir, pattern="exception.txt")) == 0){
      cat("Exception present, skipping..\n")
      next
    }
    tmpDir <- paste(mutationDir, "tmp-geonet", sep="/")
    mkdirs(tmpDir)
    setwd(tmpDir)
    file.copy(from = paste(specificGeneDir, pdbFile, sep="/"), to = tmpDir)
    system(paste(foldxExec, " --command=PositionScan --pdb=", pdbFile, " --out-pdb=true --positions=",mutation,sep=""), intern = TRUE)
    mutantAA1 <- substr(mutation, nchar(mutation), nchar(mutation))
    position <- substr(mutation, 3, nchar(mutation)-1)
    mutantAA3 <- aa1to3(mutantAA1)
    mutantPDB <- paste0(mutantAA3, position, "_", pdbFile)
    if(!file.exists(mutantPDB)){
      Sys.sleep(1)
      system(paste(foldxExec, " --command=PositionScan --pdb=", pdbFile, " --out-pdb=true --positions=",mutation,sep=""), intern = TRUE)
    }
    if(!file.exists(mutantPDB)){stop("Missing mutant pdb")}
    file.rename(mutantPDB, paste0(mutantPDB,".unfixed"))
    system(paste(pdbFixerPath, " ", paste0(mutantPDB,".unfixed"), " --output ", mutantPDB, sep=""), intern = TRUE)
    setwd(geoNetScriptsDir)
    geoNetExe_Mut <- paste0(geneName,"_", mutantAA3, "_", position, "_", geoNetExe)
    file.copy(from = geoNetExe, to = geoNetExe_Mut)
    oldArgs <- "    args = parse_args()"
    newArgs_int <- paste0("    args = parse_args([\"--querypath\", \"",tmpDir,"\", \"--filename\", \"",mutantPDB,"\", \"--chainid\", \"A\", \"--ligand\", \"P,DNA,RNA\", \"--cpu\", \"",nrOfCores,"\"])")
    ok <- replace_line_in_file(geoNetExe_Mut, oldArgs, newArgs_int)
    if(!ok){ stop("Args replace failed") }
    system(paste(pythonExe, geoNetExe_Mut, sep=" "), intern = TRUE)
    unlink(geoNetExe_Mut)
    setwd(tmpDir)
    dnaFile <- paste0(tmpDir, "/", list.files(pattern = "/*DNA-binding_result.csv"))
    rnaFile <- paste0(tmpDir, "/", list.files(pattern = "/*RNA-binding_result.csv"))
    pFile <- paste0(tmpDir, "/", list.files(pattern = "/*P-binding_result.csv"))
    combinedFile <- "CombinedGeoNetPredictions.csv"
    combine_geonet_results(dnaFile, rnaFile, pFile, c("dnaProb", "dnaBin", "rnaProb", "rnaBin", "pProb", "pBin"), combinedFile)
    if(!file.exists(combinedFile)){ stop("Combined file for mutant failed") }
    file.copy(from = combinedFile, to = mutationDir)
    if(!file.exists(paste(mutationDir, combinedFile, sep="/"))){ stop("Combined file for mutant not in final location") }
    setwd(mutationDir)
    unlink(tmpDir, recursive = TRUE)
  }
}


