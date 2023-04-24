build_peptide_df = function(y, x){
  PEPTIDE_DF = data.frame(Y = y,
                          EC = x$`Protein Accession`,
                          QValue = x$QValue,
                          DCT = x$`Decoy/Contaminant/Target`,
                          PEP = x$PEP,
                          sequence = x$`Base Sequence`)
  
  PEPTIDE_DF = PEPTIDE_DF[PEPTIDE_DF$DCT == "T",]
  
  PEPTIDE_DF
}