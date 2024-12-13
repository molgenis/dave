combine_last_columns <- function(file1, file2, file3, colNames, output_file) {
  # Read each file
  data1 <- read.csv(file1)
  data2 <- read.csv(file2)
  data3 <- read.csv(file3)
  
  # Keep only the last column of each
  last_col1 <- data1[[ncol(data1)]]
  last_col2 <- data2[[ncol(data2)]]
  last_col3 <- data3[[ncol(data3)]]
  
  # Combine the last columns into a new data frame
  combined_data <- data.frame(File1_Last_Column = last_col1,
                              File2_Last_Column = last_col2,
                              File3_Last_Column = last_col3)
  colnames(combined_data) <- colNames
  
  # Write the combined data frame to a new file
  write.csv(combined_data, file = output_file, row.names = FALSE)
}

combine_geonet_results <- function(file1, file2, file3, colNames, output_file) {
 
  # Read each file
  data1 <- read.csv(file1)
  data2 <- read.csv(file2)
  data3 <- read.csv(file3)
  
  # Keep probability and binary results
  probBin1 <- data1[,4:5]
  probBin2 <- data2[,4:5]
  probBin3 <- data3[,4:5]
  
  # Combine the last columns into a new data frame
  combined_data <- data.frame(File1_Last_Column = probBin1,
                              File2_Last_Column = probBin2,
                              File3_Last_Column = probBin3)
  colnames(combined_data) <- colNames
  
  # Write the combined data frame to a new file
  write.csv(combined_data, file = output_file, row.names = FALSE)
}
