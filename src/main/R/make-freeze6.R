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
allProtFeat <- data.frame()
for(i in seq_along(uniprotIDs))
{
  id <- uniprotIDs[i]
# id <- "P25440" # bug ones for dev/debug: P04637, P00451, P13569, P51787, P49768, Q12809, P38398, P02452
  
  # load the data and convert JSON into feature data frame
  cat(paste0("Parsing ", id, ", " , i, " of ", length(uniprotIDs), "\n"))
  inputFile <- paste0(protFeatLoc, "/", id, ".json.gz")
  protdat <- fromJSON(gzfile(inputFile))
  pf <- protdat$features
  
  # add the UniProtID
  pf$UniProtID <- id
  
  # rearrange to correct hierarchy of UniProtID -> Category -> Type
  pf <- pf %>% relocate(category, .before = type)
  pf <- pf %>% relocate(UniProtID, .before = category)
  
  # exclude rows
  pf <- pf[pf$type != "CHAIN", ] # has nice descriptions though..
  pf <- pf[pf$type != "CONFLICT", ]
  pf <- pf[pf$category != "VARIANTS", ]
  pf <- pf[pf$category != "MUTAGENESIS", ]
  
  # exclude columns
  pf$ftId <- NULL
  pf$evidences <- NULL
  pf$alternativeSequence <- NULL
  
  # override ligand or ligand part info with just ligand name
  if ("ligandPart" %in% names(pf)) {
    pf$ligand <- pf$ligandPart$name
  }else{
    pf$ligand <- pf$ligand$name
  }


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
  # move ligand over to description - implicit check on there being no ligands outside 'BINDING' type
  descr_ok_to_replace <- c("", "in other chain", "in inhibited form", "via 3-oxoalanine")
  rows <- pf$type == "BINDING" & pf$description %in% descr_ok_to_replace & !pf$ligand==""
  if(length(rows) > 0){
    pf$description[rows] <- pf$ligand[rows]
    pf$ligand[rows] <- NA
  }

  # check if 'ligand' only has NA values now and remove
  if ("ligand" %in% names(pf)) {
    if (any(!is.na(pf$ligand))) {
      stop("Error: 'ligand' column contains non-empty values.")
    } else {
      pf$ligand <- NULL
    }
  }
  
  # add to all protein features
  allProtFeat <- rbind(allProtFeat, pf)
}
allProtFeat$begin <- as.numeric(allProtFeat$begin)
allProtFeat$end <- as.numeric(allProtFeat$end)


########################################################
# Annotate freeze5 data with protein sequence features #
########################################################
frz5$seqFt <- ""
for(i in 1:nrow(frz5))
{
  i <- 1
  row <- frz5[i,]
  pos <- as.numeric(substr(row$delta_aaSeq, 3, nchar(row$delta_aaSeq)-1))
  
  allProtFeat[allProtFeat$UniProtID == "P01023"]
  
  hits <- allProtFeat[allProtFeat$UniProtID == row$UniProtID & allProtFeat$begin <= pos & allProtFeat$end >= pos]
  
}




