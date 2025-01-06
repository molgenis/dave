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


######################
# Feature extraction #
######################

# Extract GeoNet DNA/RNA/Protein binding site features from input data frames
extract_geonet_features <- function(geoNet_WT, geoNet_Mu){
  return(
    data.frame(
      delta_DNAs_cumu_prob = sum(geoNet_Mu$dnaProb) - sum(geoNet_WT$dnaProb),
      delta_DNAs_cumu_bin = sum(geoNet_Mu$dnaBin) - sum(geoNet_WT$dnaBin),
      delta_RNAs_cumu_prob = sum(geoNet_Mu$rnaProb) - sum(geoNet_WT$rnaProb),
      delta_RNAs_cumu_bin = sum(geoNet_Mu$rnaBin) - sum(geoNet_WT$rnaBin),
      delta_ProtS_cumu_prob = sum(geoNet_Mu$pProb) - sum(geoNet_WT$pProb),
      delta_ProtS_cumu_bin = sum(geoNet_Mu$pBin) - sum(geoNet_WT$pBin)
    )
  )
}

# TEMPORARY until GeoNet results are complete, bell curve values to impute for LP/LB/VUS
get_GeoNet_bellcurve_values <- function(dataGenesDir, vkglProtVarFileName)
{
  dfForImputation <- data.frame()
  
  # Iterate over selection of genes
  for(i in seq_along(succesfulGenes))
  {
    # i <- 1 # DEBUG/DEV
    geneName <- succesfulGenes[i]
    cat(paste("GeoNet imputation, working on gene:", geneName, "\n", sep=" "))
    specificGeneDir <- paste(dataGenesDir, geneName, sep="/")
    if(!length(list.files(specificGeneDir, pattern="CombinedGeoNetPredictions.csv")) == 0){
      geoNet_WT <- read.csv(paste(specificGeneDir, "CombinedGeoNetPredictions.csv", sep = "/"))
    }else{
      next
    }
    variants <- read.table(file=paste(specificGeneDir, vkglProtVarFileName, sep="/"), sep = '\t', header = TRUE, colClasses = c("character", "character", "numeric", "character", "character", "character", "character", "numeric", "character"))
    foldingResultsDir <- paste(specificGeneDir, "folding-results", sep="/")
    for(j in 1:nrow(variants))
    {
      # j <- 1 # DEBUG/DEV
      mutation <- variants$ProtChange[j]
      cat(paste("GeoNet imputation, working on ", mutation, " (gene ",geneName,", mutation ", j, " of ", nrow(variants), ")\n", sep=""))
      mutationDir <- paste(foldingResultsDir, mutation, sep="/")
      if(!length(list.files(mutationDir, pattern="CombinedGeoNetPredictions.csv")) == 0){
        geoNet_Mu <- read.csv(paste(mutationDir, "CombinedGeoNetPredictions.csv", sep="/"))
        dfForImputation <- rbind(dfForImputation, data.frame(label = variants$Classification[j], extract_geonet_features(geoNet_WT, geoNet_Mu)))
      }
    }
  }
  means <- data.frame(bell="means", aggregate(.~label, dfForImputation, mean))
  SDs <- data.frame(bell="SDs", aggregate(.~label, dfForImputation, sd))
  res <- rbind(means, SDs)
  res$bell <- as.factor(res$bell)
  res$label <- as.factor(res$label)
  return(res)
}

# TEMPORARY until GeoNet results are complete, impute for LP/LB/VUS
# # e.g. gnImp[gnImp$label=="LB" & gnImp$bell=="means",]$delta_DNAs_cumu_prob

impute_GeoNet_features <- function(gnImp, cLabel)
{
  return(
    data.frame(
      delta_DNAs_cumu_prob = rnorm(1, gnImp[gnImp$label==cLabel & gnImp$bell=="means",]$delta_DNAs_cumu_prob, gnImp[gnImp$label==cLabel & gnImp$bell=="SDs",]$delta_DNAs_cumu_prob),
      delta_DNAs_cumu_bin = rnorm(1, gnImp[gnImp$label==cLabel & gnImp$bell=="means",]$delta_DNAs_cumu_bin, gnImp[gnImp$label==cLabel & gnImp$bell=="SDs",]$delta_DNAs_cumu_bin),
      delta_RNAs_cumu_prob = rnorm(1, gnImp[gnImp$label==cLabel & gnImp$bell=="means",]$delta_RNAs_cumu_prob, gnImp[gnImp$label==cLabel & gnImp$bell=="SDs",]$delta_RNAs_cumu_prob),
      delta_RNAs_cumu_bin = rnorm(1, gnImp[gnImp$label==cLabel & gnImp$bell=="means",]$delta_RNAs_cumu_bin, gnImp[gnImp$label==cLabel & gnImp$bell=="SDs",]$delta_RNAs_cumu_bin),
      delta_ProtS_cumu_prob = rnorm(1, gnImp[gnImp$label==cLabel & gnImp$bell=="means",]$delta_ProtS_cumu_prob, gnImp[gnImp$label==cLabel & gnImp$bell=="SDs",]$delta_ProtS_cumu_prob),
      delta_ProtS_cumu_bin = rnorm(1, gnImp[gnImp$label==cLabel & gnImp$bell=="means",]$delta_ProtS_cumu_bin, gnImp[gnImp$label==cLabel & gnImp$bell=="SDs",]$delta_ProtS_cumu_bin)
    )
  )
}

# Extract GLM-Score DNA/RNA binding features from input data frames
extract_glmscore_features <- function(glmScore_WT, glmScore_Mu){
  return(
    data.frame(
      delta_DNAb_total_hydrophobic_contact_score_V2 = glmScore_Mu["V2",]$DNA - glmScore_WT["V2",]$DNA,
      delta_DNAb_Van_der_Waals_interactions_V3 = glmScore_Mu["V3",]$DNA - glmScore_WT["V3",]$DNA,
      delta_DNAb_side_chain_rotation_V4 = glmScore_Mu["V4",]$DNA - glmScore_WT["V4",]$DNA,
      delta_DNAb_hydrogen_bonding_V5 = glmScore_Mu["V5",]$DNA - glmScore_WT["V5",]$DNA,
      delta_DNAb_accessible_to_solvent_area_of_protein_V6= glmScore_Mu["V6",]$DNA - glmScore_WT["V6",]$DNA,
      delta_DNAb_accessible_to_solvent_area_of_ligand_V7 = glmScore_Mu["V7",]$DNA - glmScore_WT["V7",]$DNA,
      delta_DNAb_repulsive_interactions_V18 = glmScore_Mu["V18",]$DNA - glmScore_WT["V18",]$DNA,
      delta_DNAb_london_disperson_forces_V19 = glmScore_Mu["V19",]$DNA - glmScore_WT["V19",]$DNA,
      delta_DNAb_contact_hydrophobicity_V20 = glmScore_Mu["V20",]$DNA - glmScore_WT["V20",]$DNA,
      delta_DNAb_total_hydrophobicity_V21 = glmScore_Mu["V21",]$DNA - glmScore_WT["V21",]$DNA,
      delta_DNAb_contact_surface_tension_V22 = glmScore_Mu["V22",]$DNA - glmScore_WT["V22",]$DNA,
      delta_DNAb_total_surface_tension_V23 = glmScore_Mu["V23",]$DNA - glmScore_WT["V23",]$DNA,
      delta_DNAb_binding_affinity_pKd = glmScore_Mu["pKd",]$DNA - glmScore_WT["pKd",]$DNA,
      delta_RNAb_total_hydrophobic_contact_score_V2 = glmScore_Mu["V2",]$RNA - glmScore_WT["V2",]$RNA,
      delta_RNAb_Van_der_Waals_interactions_V3 = glmScore_Mu["V3",]$RNA - glmScore_WT["V3",]$RNA,
      delta_RNAb_side_chain_rotation_V4 = glmScore_Mu["V4",]$RNA - glmScore_WT["V4",]$RNA,
      delta_RNAb_hydrogen_bonding_V5 = glmScore_Mu["V5",]$RNA - glmScore_WT["V5",]$RNA,
      delta_RNAb_accessible_to_solvent_area_of_protein_V6 = glmScore_Mu["V6",]$RNA - glmScore_WT["V6",]$RNA,
      delta_RNAb_accessible_to_solvent_area_of_ligand_V7 = glmScore_Mu["V7",]$RNA - glmScore_WT["V7",]$RNA,
      delta_RNAb_repulsive_interactions_V18 = glmScore_Mu["V18",]$RNA - glmScore_WT["V18",]$RNA,
      delta_RNAb_london_disperson_forces_V19 = glmScore_Mu["V19",]$RNA - glmScore_WT["V19",]$RNA,
      delta_RNAb_contact_hydrophobicity_V20 = glmScore_Mu["V20",]$RNA - glmScore_WT["V20",]$RNA,
      delta_RNAb_total_hydrophobicity_V21 = glmScore_Mu["V21",]$RNA - glmScore_WT["V21",]$RNA,
      delta_RNAb_contact_surface_tension_V22 = glmScore_Mu["V22",]$RNA - glmScore_WT["V22",]$RNA,
      delta_RNAb_total_surface_tension_V23 = glmScore_Mu["V23",]$RNA - glmScore_WT["V23",]$RNA,
      delta_RNAb_binding_affinity_pKd = glmScore_Mu["pKd",]$RNA - glmScore_WT["pKd",]$RNA
      )
    )
}

# Extract P2Rank ligand pocket features from input data frames
extract_p2rank_features <- function(p2rank_WT, p2rank_Mu){
  res <- data.frame(
    delta_ligand_nr_of_predicted_pockets = dim(p2rank_Mu)[1]-dim(p2rank_WT)[1],
    delta_ligand_cumu_score = sum(p2rank_Mu$score)-sum(p2rank_WT$score),
    delta_ligand_cumu_prob = sum(p2rank_Mu$probability)-sum(p2rank_WT$probability),
    delta_ligand_cumu_sas_points = sum(p2rank_Mu$sas_points)-sum(p2rank_WT$sas_points),
    delta_ligand_cumu_surf_atoms = sum(p2rank_Mu$surf_atoms)-sum(p2rank_WT$surf_atoms),
    delta_ligand_rank1_score = p2rank_Mu$score[1]-p2rank_WT$score[1],
    delta_ligand_rank1_prob = p2rank_Mu$probability[1]-p2rank_WT$probability[1],
    delta_ligand_rank1_sas_points = p2rank_Mu$sas_points[1]-p2rank_WT$sas_points[1],
    delta_ligand_rank1_surf_atoms = p2rank_Mu$surf_atoms[1]-p2rank_WT$surf_atoms[1]
    )
  res[is.na(res)] <- 0 # if WT or Mu had 0 pockets, we get NA values. Replace with 0 here.
  return(res)
}

