from_MM_to_SE = function(path_to_peptides_psm,
                         path_to_peptides_intensities = NULL,
                         abundance_type,
                         PEP = TRUE,
                         FDR_thd = NULL) {
  
  input_type = "metamorpheus"
  variables = c("Protein Accession", "QValue", "Decoy/Contaminant/Target",
                "Base Sequence")
  if (PEP) {
    variables = c(variables, "PEP")
  }
  if (abundance_type == "intensities") {
    data_list = list(x = fread(path_to_peptides_intensities, data.table = FALSE),
                     y = NULL)
    data_list = build_intensity(data_list$x, path_to_peptides_psm,
                                variables)
  } else if (abundance_type == "psm") {
    data_list = list(x = fread(path_to_peptides_psm, select = variables,
                               data.table = FALSE),
                     y = NULL
                     )
    data_list$y = as.numeric(fread(path_to_peptides_psm,
                                   select = "PSM Count (unambiguous, <0.01 q-value)",
                                   data.table = FALSE)[, 1]
                             )
  }
  PEPTIDE_DF = do.call("build_peptide_df", data_list)
  message("We found:")
  protein_name = get_prot_from_EC(PEPTIDE_DF$EC)
  rm(data_list)
  
  PEPTIDE_DF = PEPTIDE_DF[PEPTIDE_DF$QValue <= FDR_thd,]
  
  convert_peptideDF_to_SE(PEPTIDE_DF, input_type, PEP, FDR_thd, protein_name)
}