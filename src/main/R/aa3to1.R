# Function to replace 3 letter AA codes by 1 letter
aa3to1 <- function(aa3list) {
  aa1list <- c()
  for (aa3 in aa3list) {
    aa1 <- switch(aa3, ALA={"A"}, ARG={"R"}, ASN={"N"}, ASP={"D"}, ASX={"B"}, CYS={"C"}, GLN={"Q"}, GLU={"E"}, GLX={"Z"}, GLY={"G"}, HIS={"H"}, ILE={"I"}, LEU={"L"}, LYS={"K"}, MET={"M"}, PHE={"F"}, PRO={"P"}, SER={"S"}, THR={"T"}, TRP={"W"}, TYR={"Y"}, VAL={"V"})
    aa1list <- c(aa1list, aa1)
  }
  aa1listColl <- paste(aa1list, collapse='')
  return(aa1listColl)
}

aa1to3 <- function(aa1) {
  aa3 <- switch(aa1, A={"ALA"}, R={"ARG"}, N={"ASN"}, D={"ASP"}, B={"ASX"}, C={"CYS"}, Q={"GLN"}, E={"GLU"}, Z={"GLX"}, G={"GLY"}, H={"HIS"}, I={"ILE"}, L={"LEU"}, K={"LYS"}, M={"MET"}, F={"PHE"}, P={"PRO"}, S={"SER"}, T={"THR"}, W={"TRP"}, Y={"TYR"}, V={"VAL"})
  return(aa3)
}
