# simple logistic functions
sigmoid <- function(x) { 1 / (1 + exp(-x)) }
inverseSigmoid <- function(x) { log(x / (1 - x)) }

# convert a SHAP value dataframe to SHAP probability heuristic
convert_to_SHAP_probability_heuristic <- function(df) {
  
  # Separate out the base_value (bias) from the per-feature SHAP
  base_value <- df$baseline
  shap_values <- df$shapley_values
  
  # Convert base_value to probabilities
  base_prob <- sigmoid(base_value)
  
  # Sum of all SHAP log-odds for each row
  shap_sum <- rowSums(shap_values)
  final_log_odds <- base_value + shap_sum
  final_prob <- sigmoid(final_log_odds)
  
  # Scale the per-feature SHAP by fraction of total shap_sum
  # Add a small epsilon if shap_sum is zero to avoid dividing by zero
  eps <- 1e-15
  shap_values_scaled <- shap_values / (shap_sum + eps)
  
  # Probability difference from base to final
  prob_diff <- final_prob - base_prob
  
  # Multiply scaled shap by the probability difference
  shap_values_scaled_prob <- shap_values_scaled * prob_diff
  
  # Build output dataframe
  # Rename feature columns with ".sph"
  colnames(shap_values_scaled_prob) <- paste0(colnames(shap_values), ".sph")
  result <- as.data.frame(shap_values_scaled_prob)
  
  # Add columns for base and final probabilities
  result$BaseProbability.sph  <- base_prob
  result$FinalProbability.sph <- final_prob
  
  return(result)
}