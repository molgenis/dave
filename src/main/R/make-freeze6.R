library(crunch)     # To compress results
library(R.utils)    # for 'gunzip', 'mkdirs'
library(jsonlite)   # Reading JSON into dataframes
library(dplyr)      # Data manipulation
library(tidyr)      # Data manipulation

#######################
# Adjustable settings #
#######################
rootDir <- "/Users/joeri/git/dave1" # root directory that contains README.md, data/, img/, out/, src/, etc.


########################
# Seperator characters #
########################
CONCAT_SEP <- "~" # used to combine UniProt hierarchy: category -> type -> description -> etc
VALUE_SEP <- "|" # used to combine multiple annotation values in the output data freeze

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
t(table(allProtFeat$category, allProtFeat$type))

# Sanity checks, we don't want future seperator characters in the data
rows_with_CONCAT_SEP <- rowSums(sapply(allProtFeat, grepl, pattern = CONCAT_SEP, fixed = TRUE)) > 0
if(sum(rows_with_CONCAT_SEP)>0){ stop("data contains CONCAT_SEP, print with -> allProtFeat[rows_with_CONCAT_SEP,]") }

rows_with_VALUE_SEP <- rowSums(sapply(allProtFeat, grepl, pattern = VALUE_SEP, fixed = TRUE)) > 0
if(sum(rows_with_VALUE_SEP)>0){ stop("data contains VALUE_SEP, print with -> allProtFeat[rows_with_VALUE_SEP,]") }


########################################################
# Annotate freeze5 data with protein sequence features #
########################################################
frz5$seqFt <- ""
frz5$chain <- ""
#ALLHITS <- data.frame()
for(i in 1:nrow(frz5))
{
  cat(paste0("working on ", i, " of ", nrow(frz5), "\n"))
  row <- frz5[i,]
  pos <- as.numeric(substr(row$delta_aaSeq, 3, nchar(row$delta_aaSeq)-1))
  matchID <- allProtFeat[allProtFeat$UniProtID == row$UniProtID,]
  
  # Keep CHAIN description in its own column
  getChainInfo <- matchID[matchID$type=="CHAIN",]
  if(dim(getChainInfo)[1]==0){
    getChainInfo <- getChainInfo %>% add_row(description = "")
  }
  for(k in 1:nrow(getChainInfo)) # in some cases, more than one
  {
    if(k==1){
      frz5[i,]$chain <- getChainInfo[k,]$description
    }else{
      frz5[i,]$chain <- paste0(frz5[i,]$chain, VALUE_SEP, getChainInfo[k,]$description)
    }
  }
  matchID <- matchID[matchID$type!="CHAIN",] # remove CHAIN from further processing

  # Process other protein sequence features
  matchID_DISULFID <- matchID[matchID$type=="DISULFID",]
  matchID <- matchID[matchID$type!="DISULFID",]
  hits <- matchID[pos >= matchID$begin & pos <= matchID$end,]
  # Exception: for disulphide bonds, exact match for start and end (=the locations of bonded cysteines)
  hits_DISULFID <- matchID_DISULFID[matchID_DISULFID$begin == pos & matchID_DISULFID$end == pos,]
  hits <- rbind(hits, hits_DISULFID)
  if(dim(hits)[1] == 0){ next }
  for(j in 1:nrow(hits))
  {
   #j <- 1
    #cat(paste0(". adding hit ", j, " of ", nrow(hits), "\n"))
    hit <- hits[j,]
    
    # There must ALWAYS be category and type
    if(hit$category==""){ stop("category missing!") }
    if(hit$type==""){ stop("type missing!") }
    addThis <- paste0(hit$category, VALUE_SEP, hit$category, CONCAT_SEP, hit$type)
    if(j==1){
      frz5[i,]$seqFt <- addThis
    }else{
      frz5[i,]$seqFt <- paste0(frz5[i,]$seqFt, VALUE_SEP, addThis)
    }
    
    # if we have description or extra info, add
    if(hit$description!=""){
      frz5[i,]$seqFt <- paste0(frz5[i,]$seqFt, VALUE_SEP, hit$category, CONCAT_SEP, hit$type, CONCAT_SEP, hit$description)
    }
    if(hit$extraInfo!=""){
      frz5[i,]$seqFt <- paste0(frz5[i,]$seqFt, VALUE_SEP,  hit$category, CONCAT_SEP, hit$type, CONCAT_SEP, hit$description, CONCAT_SEP, hit$extraInfo)
    }
  }
}
# Inspect results : print counts per unique feature across all seqFt values
lists <- frz5 %>%
  select(seqFt) %>%                      # keep only the one column
  filter(!is.na(seqFt)) %>%              # drop NA rows
  separate_rows(seqFt, sep = "\\|") %>%  # explode on "|"
  mutate(seqFt = trimws(seqFt)) %>%      # remove stray spaces
  filter(seqFt != "") %>%                # drop blank tokens
  count(seqFt, name = "count", sort = TRUE) %>%   # tally + sort
  rename(token = seqFt)
print(lists, n=25)
# also for chain descriptions
lists <- frz5 %>% select(chain) %>% filter(!is.na(chain)) %>% separate_rows(chain, sep = "\\|") %>% mutate(chain = trimws(chain)) %>% filter(chain != "") %>% count(chain, name = "count", sort = TRUE) %>% rename(token = chain)
print(lists, n=25)

####################
# Save as freeze 6 #
####################
resultsFreeze6 <- paste(rootDir, "data", "freeze6.csv.gz", sep="/")
write.csv.gz(frz5, resultsFreeze6, row.names = FALSE, quote = TRUE)
orig <- head(frz5)
frz5 <- NULL

# Load back in
frz6 <- read.csv(resultsFreeze6)
loaded <- head(frz6)

# Optional: assign nice row names, for instance (not part of validation)
rownames(frz6) <- paste0(frz6$gene, "/", frz6$UniProtID, ":", frz6$delta_aaSeq)

# Validate loaded against original
if (identical(orig, loaded)) {
  message("They’re identical.")
} else {
  message("They differ:")
  print(waldo::compare(orig, loaded))
}

