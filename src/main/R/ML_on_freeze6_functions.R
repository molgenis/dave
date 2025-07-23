##########################################
# Basic preparation for machine learning #
##########################################
prepForML <- function(dataFreeze)
{
  # Select only the "labeled" set: LB and LP
  dataFreeze <- dataFreeze %>% filter(ann_classificationVKGL %in% c("LB", "LP"))
  
  # Select only features relevant for training/predicting plus label (ann_classificationVKGL) and grouping (seqFt) 
  featureSelection <- c("ann_classificationVKGL", "delta_", "seqFt")
  featureRemoval   <- c("delta_aaSeq", "delta_total.energy")
  dataFreeze <- dataFreeze %>% select(matches(paste(featureSelection, collapse="|"))) %>% select(-any_of(featureRemoval))
  
  # Convert labels into booleans (LB=FALSE, LP=TRUE)
  dataFreeze <- dataFreeze %>% mutate( ann_classificationVKGL = case_when(ann_classificationVKGL == "LB" ~ FALSE, ann_classificationVKGL == "LP" ~ TRUE), ann_classificationVKGL = as.logical(ann_classificationVKGL))
  
  return(dataFreeze)
}


##############################################################
# Find sequence features with a minimum nr of both LB and LP #
##############################################################
seqFtWithMinimumNrPerLabel <- function(dataFreeze, minimum)
{
  # Get values
  df_long <- dataFreeze %>%
    mutate(value = strsplit(as.character(seqFt), "\\|")) %>%
    tidyr::unnest_longer(value) %>%
    mutate(value = trimws(value))
  # Count per classification label (boolean after prep)
  result <- df_long %>%
    group_by(value, ann_classificationVKGL) %>%
    summarise(count = n(), .groups = 'drop') %>%
    tidyr::pivot_wider(
      names_from = ann_classificationVKGL,
      values_from = count,
      names_prefix = "count_",
      values_fill = 0
    )
  # Filter by minimum, applied to each label
  result_filtered <- result %>%
    filter(count_FALSE >= minimum, count_TRUE >= minimum)
  
  return(result_filtered$value)
}


###########################################
# Slice dataframe by one sequence feature #
###########################################
sliceBySeqFt <- function(dataFreeze, seqFtSlice)
{
  filtered_df <- dataFreeze[sapply(strsplit(dataFreeze$seqFt, "\\|"), function(x) seqFtSlice %in% x), ]
  return(filtered_df)
}


####
# Train a model
####
trainModelOn <- function(dataFreeze)
{
  dataFreeze <- sfDat #debug
  featureRemoval   <- c("seqFt")
  dataFreeze <- dataFreeze %>% select(-any_of(featureRemoval))
  dataFreeze$ann_classificationVKGL <- as.factor(dataFreeze$ann_classificationVKGL)
  
  # ???
  
  model <- glm(ann_classificationVKGL ~ ., data = dataFreeze, family = binomial)
  
  rf <- randomForest(ann_classificationVKGL~., data=dataFreeze)

  
  #X <- dataFreeze %>% select(-ann_classificationVKGL)
  #y <- dataFreeze$ann_classificationVKGL
  

  #dtrain <- xgb.DMatrix(as.matrix(X), label = as.numeric(y))
  #cv <- xgb.cv(data = dtrain, nfold = 5, nrounds = 100, objective = "binary:logistic", metrics = "logloss", verbose = TRUE, early_stopping_rounds = 10)

  
}
