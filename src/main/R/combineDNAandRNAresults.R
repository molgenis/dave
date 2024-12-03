transpose_and_label <- function(dna_result_file, rna_result_file, dna_pKd, rna_pKd, output_file) {
  
  # Read in the files
  df1 <- read.table(dna_result_file, header = TRUE, sep = "\t")
  df2 <- read.table(rna_result_file, header = TRUE, sep = "\t")
  
  # Add pKD values
  df1 <- cbind(df1, pKd=dna_pKd)
  df2 <- cbind(df2, pKd=rna_pKd)
  
  # Transpose the data frames
  transposed_df1 <- t(df1)
  transposed_df2 <- t(df2)
  
  # Convert transposed matrices back to data frames
  transposed_df1 <- as.data.frame(transposed_df1)
  transposed_df2 <- as.data.frame(transposed_df2)
  
  # Rename the columns for clarity
  colnames(transposed_df1) <- "DNA"
  colnames(transposed_df2) <- "RNA"
  
  # Combine the two data frames
  combined_df <- cbind(transposed_df1, transposed_df2)
  
  # Write out the new file
  write.table(combined_df, file = output_file, sep = ",", row.names = TRUE, quote = FALSE)
  
}
