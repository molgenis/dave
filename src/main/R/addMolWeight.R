# Function to sum all combinations of a matrix diagonal and output results in a new matrix
addMolWeight <- function(collapsePDB) {
  weights <- c()
  for (i in 1:nrow(collapsePDB)) {
    resDa <- switch(collapsePDB[i,"resname"], ALA={89}, ARG={174}, ASN={132}, ASP={133}, ASX={133}, CYS={121}, GLN={146}, GLU={147}, GLX={147}, GLY={75}, HIS={155}, ILE={131}, LEU={131}, LYS={146}, MET={149}, PHE={165}, PRO={115}, SER={105}, THR={119}, TRP={204}, TYR={181}, VAL={117})
    weights <- c(weights, resDa)
  }
  weights
  collapsePDB$resDa <- weights
  return(collapsePDB)
}
