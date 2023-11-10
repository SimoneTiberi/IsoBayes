from_UserData_to_SE = function(path_to_peptides_psm,
                               abundance_type,
                               PEP = TRUE,
                               FDR_thd = NULL) {
  input_type = "other"
  if (is.data.frame(path_to_peptides_psm)) {
    PEPTIDE_DF = path_to_peptides_psm
  } else {
    PEPTIDE_DF = as.data.frame(fread(path_to_peptides_psm))
  }
  check_variables(PEPTIDE_DF)
  rm(path_to_peptides_psm)
  
  if (!is.null(PEPTIDE_DF$FDR)) {
    PEPTIDE_DF = PEPTIDE_DF[PEPTIDE_DF$FDR < FDR_thd,]
  }
  if (PEP & !("PEP" %in% colnames(PEPTIDE_DF))) {
    stop( glue("'PEP'=TRUE, but PEP not present in the data.
                      If PEP not available, please set 'PEP'=FALSE.")
          )
  }
  message("We found:")
  protein_name = get_prot_from_EC(PEPTIDE_DF$EC)
  
  convert_peptideDF_to_SE(PEPTIDE_DF, input_type, PEP, FDR_thd, protein_name)
}