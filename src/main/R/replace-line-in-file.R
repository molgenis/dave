replace_line_in_file <- function(file_name, target_string, replacement_string){
  
  # Step 1: Read the file into R
  lines <- readLines(file_name, warn = FALSE)
  
  # Step 2: Find the line matching the target string and replace it
  line_index <- grep(target_string, lines) # Find the line(s) containing the target string
  
  if (length(line_index) > 0) {
    # Replace the matching line(s) with the new string
    lines[line_index] <- replacement_string
  } else {
    warning("No line matches the specified string. No changes were made.")
    return(FALSE)
  }
  
  # Step 3: Write the modified content back to the file
  writeLines(lines, file_name)
  return(TRUE)
}
