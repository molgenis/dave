######################
# Load used packages #
######################
library(R.utils)   # for 'gunzip', 'mkdirs'
library(ggplot2)   # for plotting


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


##########################################
# Derived paths of directories and files #
##########################################
dataGenesDir <- paste(rootDir, "data", "genes", sep="/")
geneToUniprotLoc <- paste(rootDir, "data", "protein-atlas-secreted-genenames-mane-uniprot.txt", sep="/")


#################################################
# Retrieve mapping of HGNC symbol to UniProt ID #
#################################################
geneToUniprot <- read.table(file=geneToUniprotLoc, sep = '\t',header = TRUE)
head(geneToUniprot)


######################################################
# Iterate over all gene directories, repair and fold #
######################################################
setwd(dataGenesDir)
geneNames <- list.files(pattern="*", recursive=FALSE, include.dirs=TRUE)
head(geneNames)

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
  uniProtID <- geneToUniprot$UniProtKB.Swiss.Prot.ID[geneToUniprot$Gene.name==geneName]
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


##################
# Gather results #
##################
setwd(dataGenesDir)
geneNames <- list.files(pattern="*", recursive=FALSE, include.dirs=TRUE)
head(geneNames)
results <- data.frame()
for(i in seq_along(geneNames))
{
  geneName <- geneNames[i]
  cat(paste("Loading data for ", geneName, " (gene ", i, " of ", length(geneNames), ")\n", sep=""))
  specificGeneDir <- paste(dataGenesDir, geneName, sep="/")
  setwd(specificGeneDir)
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
    result$gene <- variants[j, "Gene"]
    result$protChange <- variants[j, "ProtChange"]
    result$classificationVKGL <- variants[j, "Classification"]
    results <- rbind(results, result)
  }
}

# some quick checks, replace later
a <- subset(results, classificationVKGL == "LP")
b <- subset(results, classificationVKGL == "LB")
median(a$total.energy)
median(b$total.energy)
c <- subset(results, classificationVKGL == "LP" | classificationVKGL == "LB")
ggplot(c, aes(x=total.energy, color=classificationVKGL, fill=classificationVKGL)) +
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity",binwidth=1)+
  geom_density(alpha=.2)
