build_peptide_df = function(y, x) {
  PEPTIDE_DF = data.frame(Y = y, x)
  setnames(PEPTIDE_DF, c("Protein.Accession", "Decoy.Contaminant.Target", "Base.Sequence"), c("EC", "DCT", "sequence"))
  PEPTIDE_DF = PEPTIDE_DF[PEPTIDE_DF$DCT == "T", ]

  PEPTIDE_DF
}
