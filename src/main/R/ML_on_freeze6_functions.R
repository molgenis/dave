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



######################################
# Helper functions to find constants #
######################################
is_const_num <- function(x) is.numeric(x) && (sd(x, na.rm = TRUE) == 0 || all(is.na(x)))
is_const_cat <- function(x) {
  if (!is.factor(x) && !is.character(x)) return(FALSE)
  nlevels <- length(unique(stats::na.omit(as.character(x))))
  nlevels <= 1
}


#################
# Train a model #
#################
trainModelOn <- function(dataFreeze)
{
  # Drop columns with constant values
  const_cols <- names(Filter(identity, lapply(dataFreeze, function(col) is_const_num(col) || is_const_cat(col))))
  dataFreeze <- dataFreeze[ , setdiff(names(dataFreeze), const_cols), drop = FALSE]
  
  # Remove the sequence feature for training and factorize label
  featureRemoval   <- c("seqFt")
  dataFreeze <- dataFreeze %>% select(-any_of(featureRemoval))
  dataFreeze$ann_classificationVKGL <- as.factor(dataFreeze$ann_classificationVKGL)

  # Run Bayesian generalized linear models via Stan, slow but better for small data sets
  model <- rstanarm::stan_glm(ann_classificationVKGL ~ ., data = dataFreeze, family = binomial())
  return(model)
}
