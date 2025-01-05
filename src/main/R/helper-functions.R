# Replace 3 letter AA codes by 1 letter codes
aa3to1 <- function(aa3list) {
  aa1list <- c()
  for (aa3 in aa3list) {
    aa1 <- switch(aa3, ALA={"A"}, ARG={"R"}, ASN={"N"}, ASP={"D"}, ASX={"B"}, CYS={"C"}, GLN={"Q"}, GLU={"E"}, GLX={"Z"}, GLY={"G"}, HIS={"H"}, ILE={"I"}, LEU={"L"}, LYS={"K"}, MET={"M"}, PHE={"F"}, PRO={"P"}, SER={"S"}, THR={"T"}, TRP={"W"}, TYR={"Y"}, VAL={"V"})
    aa1list <- c(aa1list, aa1)
  }
  aa1listColl <- paste(aa1list, collapse='')
  return(aa1listColl)
}

# Replace 1 letter AA codes by 3 letter codes
aa1to3 <- function(aa1) {
  aa3 <- switch(aa1, A={"ALA"}, R={"ARG"}, N={"ASN"}, D={"ASP"}, B={"ASX"}, C={"CYS"}, Q={"GLN"}, E={"GLU"}, Z={"GLX"}, G={"GLY"}, H={"HIS"}, I={"ILE"}, L={"LEU"}, K={"LYS"}, M={"MET"}, F={"PHE"}, P={"PRO"}, S={"SER"}, T={"THR"}, W={"TRP"}, Y={"TYR"}, V={"VAL"})
  return(aa3)
}

# In a target file, replace a target string with a replacement string
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

# For GLM-Score, combine the results of DNA and RNA binding (features + pKd) into one file
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

# For ScanNet output, combine the last columns of 3 result files (interface, epitope, disordered) into a one file
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

# For GeoNet, combine two columns from 3 result files (DNA, RNA, P) into one new file
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

# TEMPORARY until GeoNet results are complete, impute for LP/LB/VUS
get_imputed_GeoNet_values <- function(dataGenesDir, vkglProtVarFileName)
{
  dfForImputation <- data.frame(label=character(), delta_cumuDNAProb=numeric(), delta_cumuDNABin=numeric(), delta_cumuRNAProb=numeric(), delta_cumuRNABin=numeric(), delta_cumuPProb=numeric(), delta_cumuPBin=numeric())
  
  # Iterate over selection of genes
  for(i in seq_along(succesfulGenes))
  {
    #i <- 1 # DEBUG/DEV
    geneName <- succesfulGenes[i]
    cat(paste("GeoNet imputation, working on gene:", geneName, "\n", sep=" "))
    specificGeneDir <- paste(dataGenesDir, geneName, sep="/")
    if(!length(list.files(specificGeneDir, pattern="CombinedGeoNetPredictions.csv")) == 0){
      wt_gnRes <- read.csv(paste(specificGeneDir, "CombinedGeoNetPredictions.csv", sep = "/"))
    }else{
      next
    }
    
    wt_cumuDNAProb <- sum(wt_gnRes$dnaProb)
    wt_cumuDNABin <- sum(wt_gnRes$dnaBin)
    wt_cumuRNAProb <- sum(wt_gnRes$rnaProb)
    wt_cumuRNABin <- sum(wt_gnRes$rnaBin)
    wt_cumuPProb <- sum(wt_gnRes$pProb)
    wt_cumuPBin <- sum(wt_gnRes$pBin)
    
    variants <- read.table(file=paste(specificGeneDir, vkglProtVarFileName, sep="/"), sep = '\t', header = TRUE, colClasses = c("character", "character", "numeric", "character", "character", "character", "character", "numeric", "character"))
    foldingResultsDir <- paste(specificGeneDir, "folding-results", sep="/")
    
    for(j in 1:nrow(variants))
    {
      #j <- 1 # DEBUG/DEV
      mutation <- variants$ProtChange[j]
      cat(paste("GeoNet imputation, working on ", mutation, " (gene ",geneName,", mutation ", j, " of ", nrow(variants), ")\n", sep=""))
      mutationDir <- paste(foldingResultsDir, mutation, sep="/")
      if(!length(list.files(mutationDir, pattern="CombinedGeoNetPredictions.csv")) == 0){
        mut_gnRes <- read.csv(paste(mutationDir, "CombinedGeoNetPredictions.csv", sep="/"))
        mut_cumuDNAProb <- sum(mut_gnRes$dnaProb)
        mut_cumuDNABin <- sum(mut_gnRes$dnaBin)
        mut_cumuRNAProb <- sum(mut_gnRes$rnaProb)
        mut_cumuRNABin <- sum(mut_gnRes$rnaBin)
        mut_cumuPProb <- sum(mut_gnRes$pProb)
        mut_cumuPBin <- sum(mut_gnRes$pBin)
        
        dfForImputation <- rbind(dfForImputation, data.frame(label = variants$Classification[j],
                                                             delta_cumuDNAProb = mut_cumuDNAProb - wt_cumuDNAProb,
                                                             delta_cumuDNABin = mut_cumuDNABin - wt_cumuDNABin,
                                                             delta_cumuRNAProb = mut_cumuRNAProb - wt_cumuRNAProb,
                                                             delta_cumuRNABin = mut_cumuRNABin - wt_cumuRNABin,
                                                             delta_cumuPProb = mut_cumuPProb - wt_cumuPProb,
                                                             delta_cumuPBin = mut_cumuPBin - wt_cumuPBin))
      }
    }
  }
  agg <- aggregate(.~label, dfForImputation, mean)
  return(agg)
}

# Extract GLM-Score DNA/RNA binding features from input data frames
extract_glmscore_features <- function(glmScore_WT, glmScore_Mu){

  featDF <- data.frame(delta_DNA_total_hydrophobic_contact_score_V2 = numeric(),
                       delta_DNA_Van_der_Waals_interactions_V3 = numeric(),
                       delta_DNA_side_chain_rotation_V4 = numeric(),
                       delta_DNA_hydrogen_bonding_V5 = numeric(),
                       delta_DNA_accessible_to_solvent_area_of_protein_V6 = numeric(),
                       delta_DNA_accessible_to_solvent_area_of_ligand_V7 = numeric(),
                       delta_DNA_repulsive_interactions_V18 = numeric(),
                       delta_DNA_london_disperson_forces_V19 = numeric(),
                       delta_DNA_contact_hydrophobicity_V20 = numeric(),
                       delta_DNA_total_hydrophobicity_V21 = numeric(),
                       delta_DNA_contact_surface_tension_V22 = numeric(),
                       delta_DNA_total_surface_tension_V23 = numeric(),
                       delta_DNA_binding_affinity_pKd = numeric(),
                       delta_RNA_total_hydrophobic_contact_score_V2 = numeric(),
                       delta_RNA_Van_der_Waals_interactions_V3 = numeric(),
                       delta_RNA_side_chain_rotation_V4 = numeric(),
                       delta_RNA_hydrogen_bonding_V5 = numeric(),
                       delta_RNA_accessible_to_solvent_area_of_protein_V6 = numeric(),
                       delta_RNA_accessible_to_solvent_area_of_ligand_V7 = numeric(),
                       delta_RNA_repulsive_interactions_V18 = numeric(),
                       delta_RNA_london_disperson_forces_V19 = numeric(),
                       delta_RNA_contact_hydrophobicity_V20 = numeric(),
                       delta_RNA_total_hydrophobicity_V21 = numeric(),
                       delta_RNA_contact_surface_tension_V22 = numeric(),
                       delta_RNA_total_surface_tension_V23 = numeric(),
                       delta_RNA_binding_affinity_pKd = numeric())

  return(
    data.frame(
      delta_DNA_total_hydrophobic_contact_score_V2 = glmScore_Mu["V2",]$DNA - glmScore_WT["V2",]$DNA,
      delta_DNA_Van_der_Waals_interactions_V3 = glmScore_Mu["V3",]$DNA - glmScore_WT["V3",]$DNA
      )
  )
}

# Extract P2Rank ligand pocket features from input data frames
extract_p2rank_features <- function(p2rank_WT, p2rank_Mu){
  
  
  featDF <- data.frame(delta_nr_of_predicted_pockets = numeric(),
                       delta_cumu_score = numeric(),
                       delta_cumu_prob = numeric(),
                       delta_cumu_sas_points = numeric(),
                       delta_cumu_surf_atoms = numeric(),
                       delta_rank1_score = numeric(),
                       delta_rank1_prob = numeric(),
                       delta_rank1_sas_points = numeric(),
                       delta_rank1_surf_atoms = numeric()
                       )
  
  
  
}

