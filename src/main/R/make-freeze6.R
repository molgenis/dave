library(crunch)     # To compress results
library(R.utils)    # for 'gunzip', 'mkdirs'
library(jsonlite)   # Reading JSON into dataframes


#######################
# Adjustable settings #
#######################
vkglProtVarFileName <- "VKGL_apr2024_protForFolding.tsv" # loading processed VKGL protein variants from: /data/genes/{gene}/{vkglProtVarFile}
rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability" # root directory that contains README.md, data/, img/, out/, src/, etc.


################################
# Get uniprot IDs from freeze5 #
################################
frz5loc <- paste(rootDir, "data", "freeze5.csv.gz", sep="/")
frz5 <- read.csv(frz5loc)
uniprotIDs <- unique(frz5$UniProtID)


##########################################
# Download and protein sequence features #
##########################################
outputDir <- paste(rootDir, "data", "protfeat", sep="/")
mkdirs(outputDir)
for(id in uniprotIDs){
  # id <- "P01023" # for debugging
  outputFile <- paste0(protFeatLoc, "/", id, ".json.gz")
  if(file.exists(outputFile))
  {
    cat(paste0(outputFile," exists, skipping...\n"))
    next
  }
  url  <- paste0("https://www.ebi.ac.uk/proteins/api/features/", id)
  txt  <- readLines(url, encoding = "UTF-8", warn = FALSE)
  write(txt, outputFile)
  gzip(outputFile)
  Sys.sleep(1) # be nice citizens, take it slow
}


####################################################
# Read protein sequence features and make freeze 6 #
####################################################
for(id in uniprotIDs){
  id <- "P04637" # for debugging
  inputFile <- paste0(protFeatLoc, "/", id, ".json.gz")
  protdat <- fromJSON(gzfile(inputFile))
  protfeat <- protdat$features
  # exclude:
  # CATEGORY == VARIANTS, MUTAGENESIS, 
  protfeat <- protfeat[protfeat$category != "VARIANTS", ]
  protfeat <- protfeat[protfeat$category != "MUTAGENESIS", ]
  # TODO !!
  # filter more??
  # add gene, then concat to big dataframe?
}
# use big dataframe to annotate freeze5 and write out freeze 6?


