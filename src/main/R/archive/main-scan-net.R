######################
# Load used packages #
######################
library(R.utils)    # for 'gunzip', 'mkdirs'
library(reticulate) # for Python venv --> see InstallingScanNetPythonVenv.txt

#######################
# Adjustable settings #
#######################
# Example for Unix
rootDir <- "/Users/joeri/git/dave1" # root directory that contains README.md, data/, img/, out/, src/, etc.
vkglProtVarFileName <- "VKGL_apr2024_protForFolding.tsv" # loading processed VKGL protein variants from: /data/genes/{gene}/{vkglProtVarFile}
foldxExec <- "/Applications/FoldX/5/foldx5MacStd/foldx_20241231" # exact path to the FoldX executable
scanNetDir <- "/Users/joeri/git/ScanNet"
pythonVenv <- "/Users/joeri/git/ScanNet/venv"
scanNetExe <- "predict_bindingsites.py"
# Example for Windows
rootDir <- "D:/github/dave1"
vkglProtVarFileName <- "VKGL_apr2024_protForFolding.tsv"
foldxExec <- "C:/\"Program Files\"/FoldX/foldx_20241231.exe"
scanNetDir <- "D:/github/ScanNet"
pythonVenv <- "D:\\github\\ScanNet\\venv"
scanNetExe <- "predict_bindingsites.py"

# Assuming data was produced by main.R and enriched by CombineWithAlhaMissense.R
# Load the freeze4 data and assign meaningful row names
freeze <- paste(rootDir, "data", "freeze4.csv.gz", sep="/")
results <- read.csv(freeze)
rownames(results) <- paste0(results$gene, "/", results$UniProtID, ":", results$delta_aaSeq)
dataGenesDir <- paste(rootDir, "data", "genes", sep="/")
source(paste(rootDir, "src", "main", "R", "aa3to1.R", sep="/"))
source(paste(rootDir, "src", "main", "R", "replace-line-in-file.R", sep="/"))
source(paste(rootDir, "src", "main", "R", "combine-last-columns.R", sep="/"))
use_virtualenv(pythonVenv)

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
  if(!length(list.files(specificGeneDir, pattern="CombinedScanNetPredictions.csv")) == 0){
    cat("Wild-type CombinedScanNetPredictions.csv file already present, skipping..\n")
    
  }else{
    
    tmpDir <- paste(specificGeneDir, "tmp", sep="/")
    mkdirs(tmpDir)
    setwd(tmpDir)
    file.copy(from = paste(specificGeneDir, pdbFile, sep="/"), to = tmpDir)
    
    setwd(scanNetDir)
    
    scanNetExe_WT_int <- paste0(geneName,"_WT_int_", scanNetExe)
    scanNetExe_WT_epi <- paste0(geneName,"_WT_epi_", scanNetExe)
    scanNetExe_WT_idp <- paste0(geneName,"_WT_idp_", scanNetExe)
    
    file.copy(from = scanNetExe, to = scanNetExe_WT_int)
    file.copy(from = scanNetExe, to = scanNetExe_WT_epi)
    file.copy(from = scanNetExe, to = scanNetExe_WT_idp)
    
    oldArgs <- "    args = parser.parse_args()"
    newArgs_int <- paste0("    args = parser.parse_args([\"",paste(tmpDir, pdbFile, sep="/"),"\", \"--noMSA\", \"--predictions_folder\", \"",tmpDir,"\", \"--mode\", \"interface\"])")
    newArgs_epi <- paste0("    args = parser.parse_args([\"",paste(tmpDir, pdbFile, sep="/"),"\", \"--noMSA\", \"--predictions_folder\", \"",tmpDir,"\", \"--mode\", \"epitope\"])")
    newArgs_idp <- paste0("    args = parser.parse_args([\"",paste(tmpDir, pdbFile, sep="/"),"\", \"--noMSA\", \"--predictions_folder\", \"",tmpDir,"\", \"--mode\", \"idp\"])")
    
    ok1 <- replace_line_in_file(scanNetExe_WT_int, oldArgs, newArgs_int)
    ok2 <- replace_line_in_file(scanNetExe_WT_epi, oldArgs, newArgs_epi)
    ok3 <- replace_line_in_file(scanNetExe_WT_idp, oldArgs, newArgs_idp)
    
    if(!(ok1 && ok2 && ok3)){
      stop("Args replace failed")
    }
    
    py_run_file(scanNetExe_WT_int)
    py_run_file(scanNetExe_WT_epi)
    py_run_file(scanNetExe_WT_idp)
    
    unlink(scanNetExe_WT_int)
    unlink(scanNetExe_WT_epi)
    unlink(scanNetExe_WT_idp)
    
    setwd(tmpDir)
    pdbNoExt <- substr(pdbFile, 1, nchar(pdbFile)-4)
    intFile <- paste0(tmpDir, "/", pdbNoExt, "_single_ScanNet_interface_noMSA/", "predictions_",pdbNoExt,".csv")
    epiFile <- paste0(tmpDir, "/", pdbNoExt, "_single_ScanNet_epitope_noMSA/", "predictions_",pdbNoExt,".csv")
    idpFile <- paste0(tmpDir, "/", pdbNoExt, "_single_ScanNet_idp_noMSA/", "predictions_",pdbNoExt,".csv")
    combinedFile <- "CombinedScanNetPredictions.csv"
    combine_last_columns(intFile, epiFile, idpFile, c("int", "epi", "idp"), combinedFile)
    
    if(!file.exists(combinedFile)){
      stop("Combined file for WT failed")
    }
    
    file.copy(from = combinedFile, to = specificGeneDir)
    
    if(!file.exists(paste(specificGeneDir, combinedFile, sep="/"))){
      stop("Combined file for WT not in final location")
    }
    
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
    if(!length(list.files(mutationDir, pattern="CombinedScanNetPredictions.csv")) == 0){
      cat("Mutant CombinedScanNetPredictions.csv file already present, skipping..\n")
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
    
    if(!file.exists(mutantPDB)){
      Sys.sleep(1)
      system(paste(foldxExec, " --command=PositionScan --pdb=", pdbFile, " --out-pdb=true --positions=",mutation,sep=""), intern = TRUE)
    }

    setwd(scanNetDir)
    
    scanNetExe_mu_int <- paste0(geneName,"_",mutantAA3, position, "_int_", scanNetExe)
    scanNetExe_mu_epi <- paste0(geneName,"_",mutantAA3, position, "_epi_", scanNetExe)
    scanNetExe_mu_idp <- paste0(geneName,"_",mutantAA3, position, "_idp_", scanNetExe)
    
    file.copy(from = scanNetExe, to = scanNetExe_mu_int)
    file.copy(from = scanNetExe, to = scanNetExe_mu_epi)
    file.copy(from = scanNetExe, to = scanNetExe_mu_idp)
    
    oldArgs <- "    args = parser.parse_args()"
    newArgs_int <- paste0("    args = parser.parse_args([\"",paste(tmpDir, mutantPDB, sep="/"),"\", \"--noMSA\", \"--predictions_folder\", \"",tmpDir,"\", \"--mode\", \"interface\"])")
    newArgs_epi <- paste0("    args = parser.parse_args([\"",paste(tmpDir, mutantPDB, sep="/"),"\", \"--noMSA\", \"--predictions_folder\", \"",tmpDir,"\", \"--mode\", \"epitope\"])")
    newArgs_idp <- paste0("    args = parser.parse_args([\"",paste(tmpDir, mutantPDB, sep="/"),"\", \"--noMSA\", \"--predictions_folder\", \"",tmpDir,"\", \"--mode\", \"idp\"])")
    
    ok1 <- replace_line_in_file(scanNetExe_mu_int, oldArgs, newArgs_int)
    ok2 <- replace_line_in_file(scanNetExe_mu_epi, oldArgs, newArgs_epi)
    ok3 <- replace_line_in_file(scanNetExe_mu_idp, oldArgs, newArgs_idp)
    
    if(!(ok1 && ok2 && ok3)){
      stop("Args replace failed")
    }
    
    py_run_file(scanNetExe_mu_int)
    py_run_file(scanNetExe_mu_epi)
    py_run_file(scanNetExe_mu_idp)
    
    unlink(scanNetExe_mu_int)
    unlink(scanNetExe_mu_epi)
    unlink(scanNetExe_mu_idp)
    
    setwd(tmpDir)
    pdbNoExt <- substr(mutantPDB, 1, nchar(mutantPDB)-4)
    intFile <- paste0(tmpDir, "/", pdbNoExt, "_single_ScanNet_interface_noMSA/", "predictions_",pdbNoExt,".csv")
    epiFile <- paste0(tmpDir, "/", pdbNoExt, "_single_ScanNet_epitope_noMSA/", "predictions_",pdbNoExt,".csv")
    idpFile <- paste0(tmpDir, "/", pdbNoExt, "_single_ScanNet_idp_noMSA/", "predictions_",pdbNoExt,".csv")
    combinedFile <- "CombinedScanNetPredictions.csv"
    combine_last_columns(intFile, epiFile, idpFile, c("int", "epi", "idp"), combinedFile)
    
    if(!file.exists(combinedFile)){
      stop("Combined file for mutation failed")
    }
    
    file.copy(from = combinedFile, to = mutationDir)
    
    if(!file.exists(paste(mutationDir, combinedFile, sep="/"))){
      stop("Combined file for mutation not in final location")
    }

    setwd(mutationDir)
    unlink(tmpDir, recursive = TRUE)
  }
}


