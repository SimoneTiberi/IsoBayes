# function to load results, with options for:
# - PEP and no PEP inference;
# - PSM or intensities (1 or 2 file paths);
# - percolator or metamorpheus input.

# general todo:
# remove "library" and put all functions from other packages in roxygen_tags.R
#' @export
load_data = function(path_to_peptides_psm,
                      path_to_peptides_intensities = NULL,
                      tpm_path = NULL,
                      path_fasta = NULL,
                      input_type = "metamorpheus",
                      abundance_type = "psm",
                      PEP = FALSE,
                      FDR_thd = 0.01 # ignored if input_type = openMS
) {
  if (input_type == "openMS"){
    print("input_type = 'openMS'. 'FDR_thd' is ignored and abundance_type = 'psm'.")
  }
  variables = c("Protein Accession", "QValue", "Decoy/Contaminant/Target", "Base Sequence")
  if (input_type == "metamorpheus") {
    if(PEP){
      variables = c(variables, "PEP")
    }
    if (abundance_type == "intensities") {
      data_list = list(x = fread(path_to_peptides_intensities, data.table = FALSE), y = NULL)
      data_list = build_intensity(data_list$x, path_to_peptides_psm, variables)
    } else if (abundance_type == "psm") {
      data_list = list(x = fread(path_to_peptides_psm, select = variables, data.table = FALSE), y = NULL)
      data_list$y = as.numeric(fread(path_to_peptides_psm, select = "PSM Count (unambiguous, <0.01 q-value)", data.table = FALSE)[, 1])
    } else {
      stop("Invalid abundance_type Choose one of 'psm' or 'intensities'.")
    }
    PEPTIDE_DF = do.call("build_peptide_df", data_list)
    rm(data_list)
    protein_df_args = list(protein_name = get_prot_from_EC(PEPTIDE_DF$EC))
  } else if (input_type == "openMS"){
    PEPTIDE_DF = get_peptides_from_idXML(path_to_peptides_psm, PEP)
    PROTEIN_DF_openMS = get_proteins_from_idXML(path_to_peptides_psm)
    protein_name_openMS = get_prot_from_EC(PEPTIDE_DF$EC)
    protein_df_args = list(
      protein_name = PROTEIN_DF_openMS$isoform[match(protein_name_openMS, PROTEIN_DF_openMS$id)],
      id_openMS = PROTEIN_DF_openMS$id[match(protein_name_openMS, PROTEIN_DF_openMS$id)]
    )
  } else {
    stop("Invalid input_type. Choose one of 'metamorpheus' or 'openMS'.")
  }
  if (!is.null(tpm_path)) {
    protein_df_args$TPM = load_tpm(protein_df_args$protein_name, tpm_path)
  } 
  if (!is.null(path_fasta)) {
    protein_df_args$protein_length = get_protein_length(protein_df_args$protein_name, path_fasta)
  } else {
    protein_df_args$protein_length = rep(1, length(protein_df_args$protein_name))
  }
  
  PROTEIN_DF = do.call("data.frame", protein_df_args)
  rm(protein_df_args)

  if (input_type == "metamorpheus") {
    PEPTIDE_DF = PEPTIDE_DF[PEPTIDE_DF$QValue <= FDR_thd, ]
  }
  
  PEPTIDE_DF = collapse_pept_w_equal_EC(PEPTIDE_DF, PEP)
  PEPTIDE_DF = PEPTIDE_DF[PEPTIDE_DF$Y > 0, ]

  print("Total number of proteins we MAY actually detect:")
  protein_name_to_keep = get_prot_from_EC(PEPTIDE_DF$EC)
  if (input_type == "metamorpheus") {
    PROTEIN_DF = PROTEIN_DF[PROTEIN_DF$protein_name %in% protein_name_to_keep, ]
    UNIQUE_PEPT_ABUNDANCE = unique_protein_abundance(PEPTIDE_DF$Y, PEPTIDE_DF$EC, PROTEIN_DF$protein_name)
  } else {
    PROTEIN_DF = PROTEIN_DF[PROTEIN_DF$id_openMS %in% protein_name_to_keep, ]
    UNIQUE_PEPT_ABUNDANCE = unique_protein_abundance(PEPTIDE_DF$Y, PEPTIDE_DF$EC, PROTEIN_DF$id_openMS)
  }
  PROTEIN_DF$Y_unique = UNIQUE_PEPT_ABUNDANCE$Y_unique

  if (PEP) {
    PEPTIDE_DF_unique = PEPTIDE_DF[UNIQUE_PEPT_ABUNDANCE$sel_unique, ]
    PEPTIDE_DF_unique$EC_numeric = UNIQUE_PEPT_ABUNDANCE$EC_numeric[UNIQUE_PEPT_ABUNDANCE$sel_unique]
    PEPTIDE_DF_unique$EC = NULL
  } else {
    PEPTIDE_DF_unique = NULL
  }

  PEPTIDE_DF = PEPTIDE_DF[!UNIQUE_PEPT_ABUNDANCE$sel_unique, ] # keep multi-mapping peptides
  PEPTIDE_DF$EC_numeric = UNIQUE_PEPT_ABUNDANCE$EC_numeric[!UNIQUE_PEPT_ABUNDANCE$sel_unique]
  PEPTIDE_DF$EC = NULL

  print(glue("Number of multi-mapping peptides: {nrow(PEPTIDE_DF)}"))
  list(PEPTIDE_DF = PEPTIDE_DF, PEPTIDE_DF_unique = PEPTIDE_DF_unique, PROTEIN_DF = PROTEIN_DF, PEP = PEP)
}
