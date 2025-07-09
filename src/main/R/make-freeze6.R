library(crunch)     # To compress results
library(R.utils)    # for 'gunzip', 'mkdirs'
library(jsonlite)   # Reading JSON into dataframes
library(dplyr)


#######################
# Adjustable settings #
#######################
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
  outputFile <- paste0(outputDir, "/", id, ".json.gz")
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
allProtFeat <- data.frame()
for(i in seq_along(uniprotIDs))
{
  id <- uniprotIDs[i]
# id <- "P01023" # bug ones for dev/debug: P04637, P00451, P13569, P51787, P49768, Q12809, P38398, P02452 # P25440 with ligand parts
  
  # load the data and convert JSON into feature data frame
  cat(paste0("Parsing ", id, ", " , i, " of ", length(uniprotIDs), "\n"))
  inputFile <- paste0(outputDir, "/", id, ".json.gz")
  protdat <- fromJSON(gzfile(inputFile))
  pf <- protdat$features
  
  # add the UniProtID
  pf$UniProtID <- id
  
  # rearrange to correct hierarchy of UniProtID -> Category -> Type
  pf <- pf %>% relocate(category, .before = type)
  pf <- pf %>% relocate(UniProtID, .before = category)
  
  # exclude rows
  pf <- pf[pf$type != "CONFLICT", ]
  pf <- pf[pf$category != "VARIANTS", ]
  pf <- pf[pf$category != "MUTAGENESIS", ]
  
  # exclude columns
  pf$ftId <- NULL
  pf$evidences <- NULL
  pf$alternativeSequence <- NULL
  pf$molecule <- NULL
  
  # flatten and cleanup of description values:
  # simplify values containing 'cleavage' or 'interaction'
  pf$description[pf$category == "DOMAINS_AND_SITES" & grepl("cleavage", pf$description, ignore.case = TRUE)] <- "Cleavage"
  pf$description[pf$category == "DOMAINS_AND_SITES" & grepl("interaction", pf$description, ignore.case = TRUE)] <- "Interaction"
  # remove data after semicolon as it makes it unnecessarily specific
  has_semicolon <- grepl(";", pf$description, fixed = TRUE)
  pf$description[has_semicolon] <- sub(";.*$", "", pf$description[has_semicolon])
  # remove trailing numbers that indicate multiple effects/bindings of the same type
  has_trailing_num <- grepl("\\s\\d+$", pf$description)
  pf$description[has_trailing_num] <- sub("^(.*\\S)\\s+\\d+$", "\\1",pf$description[has_trailing_num], perl = TRUE)
  # override ligand or ligand part info with just their name, ignored by R if these do not exist
  pf$ligand <- pf$ligand$name
  pf$ligandPart <- pf$ligandPart$name
  
  # add ligand/ligandPart columns if missing to make concat easier
  if (!"ligand" %in% names(pf)) {
    pf$ligand <- ""
  }
  if (!"ligandPart" %in% names(pf)) {
    pf$ligandPart <- ""
  }
  # replace all NA with "" to make pasting strings easier
  pf[is.na(pf)] <- ""
  
  # for rows with type BINDING that have a description, add it to ligandPart
  # e.g. ligand 'Fe' + description 'via 3-oxoalanine' -> 'Fe via 3-oxoalanine'
  # e.g. '' + 'via 3-oxoalanine' -> 'via 3-oxoalanine'
  # e.g. 'Fe' + '' -> 'Fe'
  rows <- pf$type == "BINDING" & !pf$description=="" & !pf$ligandPart==""
  if(any(rows)){
    pf$ligandPart[rows] <- paste0(pf$ligandPart[rows], " ", pf$description[rows])
  }
  rows <- pf$type == "BINDING" & !pf$description=="" & pf$ligandPart==""
  if(any(rows)){
    pf$ligandPart[rows] <- paste0(pf$description[rows])
  }
  # then, move ligand to description
  rows <- pf$type == "BINDING" & !pf$ligand==""
  if(any(rows)){
    pf$description[rows] <- pf$ligand[rows]
  }
  # check if data other than BINDING had ligands (should NOT be the case)
  rows <- pf$type != "BINDING" & !pf$ligand==""
  if(any(rows)){
    stop("ligand has a value but type is not BINDING")
  }
  # rename ligandPart to extraInfo, remove ligand and ligandPart
  pf$ligand <- NULL
  pf$extraInfo <- pf$ligandPart
  pf$ligandPart <- NULL
  
  # add to all protein features
  allProtFeat <- rbind(allProtFeat, pf)
}
# Post-process:
# Replace missing or uncertain coordinate notation '~' with ''
allProtFeat <- allProtFeat %>% mutate(across(where(~ is.character(.x) | is.factor(.x)), ~ gsub("~", "", as.character(.x), fixed = TRUE)))
# Then remove any rows missing a begin or end (or both..)
allProtFeat <- allProtFeat[!(allProtFeat$begin == "" | allProtFeat$end == ""),]
# Then convert begin/end to numeric
allProtFeat$begin <- as.numeric(allProtFeat$begin)
allProtFeat$end <- as.numeric(allProtFeat$end)
# Make nice overview table of major effect types
table(allProtFeat$category, allProtFeat$type)


########################################################
# Annotate freeze5 data with protein sequence features #
########################################################
frz5$seqFt <- ""
ALLHITS <- data.frame()
for(i in 1:nrow(frz5))
{
# i <- 4
  cat(paste0("working on ", i, " of ", nrow(frz5), "\n"))
  row <- frz5[i,]
  pos <- as.numeric(substr(row$delta_aaSeq, 3, nchar(row$delta_aaSeq)-1))
  pos
  matchID <- allProtFeat[allProtFeat$UniProtID == row$UniProtID,]
  matchID <- matchID[matchID$type!="CHAIN",] # skip CHAIN for now
  matchID_DISULFID <- matchID[matchID$type=="DISULFID",]
  matchID <- matchID[matchID$type!="DISULFID",]

  hits <- matchID[pos >= matchID$begin & pos <= matchID$end,]
  # Exception: for disulphide bonds, exact match for start and end (=the locations of bonded cysteines)
  hits_DISULFID <- matchID_DISULFID[matchID_DISULFID$begin == pos & matchID_DISULFID$end == pos,]
  hits <- rbind(hits, hits_DISULFID)
  ALLHITS <- rbind(ALLHITS, hits)
}




