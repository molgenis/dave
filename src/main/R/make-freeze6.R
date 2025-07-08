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
  id <- "P04637" # bug ones for dev/debug: P04637, P00451, P13569, P51787, P49768, Q12809, P38398, P02452
  inputFile <- paste0(protFeatLoc, "/", id, ".json.gz")
  protdat <- fromJSON(gzfile(inputFile))
  protfeat <- protdat$features
  
  # exclude rows
  protfeat <- protfeat[protfeat$type != "CHAIN", ] # has nice descriptions though..
  protfeat <- protfeat[protfeat$type != "CONFLICT", ]
  protfeat <- protfeat[protfeat$category != "VARIANTS", ]
  protfeat <- protfeat[protfeat$category != "MUTAGENESIS", ]
  
  # exclude columns
  protfeat$ftId <- NULL
  protfeat$evidences <- NULL
  protfeat$alternativeSequence <- NULL
  
  # override ligand info with just ligand name
  protfeat$ligand <- protfeat$ligand$name

  # flatten and cleanup of description values:
  # simplify values containing 'cleavage' or 'interaction'
  protfeat$description[protfeat$category == "DOMAINS_AND_SITES" & grepl("cleavage", protfeat$description, ignore.case = TRUE)] <- "Cleavage"
  protfeat$description[protfeat$category == "DOMAINS_AND_SITES" & grepl("interaction", protfeat$description, ignore.case = TRUE)] <- "Interaction"
  # remove data after semicolon as it makes it unnecessarily specific
  has_semicolon <- grepl(";", protfeat$description, fixed = TRUE)
  protfeat$description[has_semicolon] <- sub(";.*$", "", protfeat$description[has_semicolon])
  # remove trailing numbers that indicate multiple effects/bindings of the same type
  has_trailing_num <- grepl("\\s\\d+$", protfeat$description)
  protfeat$description[has_trailing_num] <- sub("^(.*\\S)\\s+\\d+$", "\\1",protfeat$description[has_trailing_num], perl = TRUE)
  # move ligand over to description
  rows <- protfeat$type == "BINDING" & protfeat$description=="" & !protfeat$ligand==""
  protfeat$description[rows] <- protfeat$ligand[rows]  
  
  # TODO !!
  # filter more??
  # add gene, then concat to big dataframe?
}
# use big dataframe to annotate freeze5 and write out freeze 6?


