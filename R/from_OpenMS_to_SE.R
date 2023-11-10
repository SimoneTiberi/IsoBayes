from_OpenMS_to_SE = function(path_to_peptides_psm,
                             PEP = TRUE,
                             FDR_thd = NULL) {
  input_type = "openMS"
  file_idXML = fread(path_to_peptides_psm, sep = NULL, header = FALSE)
  PEPTIDE_DF = get_peptides_from_idXML(file_idXML, PEP, FDR_thd)
  message("We found:")
  PROTEIN_DF_openMS = get_proteins_from_idXML(file_idXML)
  
  protein_name_openMS = get_prot_from_EC(PEPTIDE_DF$EC)
  protein_name = PROTEIN_DF_openMS$isoform[match(protein_name_openMS,
                                                 PROTEIN_DF_openMS$id)]
  id_openMS = PROTEIN_DF_openMS$id[match(protein_name_openMS,
                                         PROTEIN_DF_openMS$id)]
  
  convert_peptideDF_to_SE(PEPTIDE_DF, input_type, PEP, FDR_thd, protein_name,
                          id_openMS)
}