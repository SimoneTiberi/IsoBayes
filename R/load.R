#' Load data
#'
#' \code{load_data} read and process input data (peptides counts from metamorpheus or OpenMS, peptides intensities, tpm)
#'
#' @param path_to_peptides_psm a character string indicating the path to psmtsv file from metamorpheus tool or the idXML file from OpenMS toolkit.
#' @param path_to_peptides_intensities a character string indicating the path to psmtsv file from metamorpheus with intensity values.
#' @param tpm_path a character string indicating the path to tpm file.
#' @param input_type a character string indicating the tool that outputs the peptides file: 'metamorpheus' or 'openMS'.
#' @param abundance_type a character string indicating the type of input: 'psm' or 'intensities'.
#' @param PEP boolean value indicating if Peptite Error Probability shoud be used.
#' @param FDR_thd numeric value indicating the False Discovery Rate filter.
#'
#' @return A list with a peptide dataframe, a peptide dataframe of unique peptides (only for PEP=TRUE) and a protein dataframe.
#' @examples
#' @author name
#'
#' @seealso links
#'
#' @export
load_data = function(path_to_peptides_psm,
                     path_to_peptides_intensities = "",
                     tpm_path = "",
                     input_type = "metamorpheus",
                     abundance_type = "psm",
                     PEP = FALSE,
                     FDR_thd = 0.01 # ignored if input_type = openMS
) {
  input_check(path_to_peptides_psm, path_to_peptides_intensities, tpm_path, input_type, abundance_type, PEP, FDR_thd)
  
  if (input_type == "openMS") {
    abundance_type = "psm"
    print("input_type = 'openMS'. 'FDR_thd' is ignored and abundance_type = 'psm'.")
  }
  variables = c("Protein Accession", "QValue", "Decoy/Contaminant/Target", "Base Sequence")
  if (input_type == "metamorpheus") {
    if (PEP) {
      variables = c(variables, "PEP")
    }
    if (abundance_type == "intensities") {
      data_list = list(x = fread(path_to_peptides_intensities, data.table = FALSE), y = NULL)
      data_list = build_intensity(data_list$x, path_to_peptides_psm, variables)
    } else if (abundance_type == "psm") {
      data_list = list(x = fread(path_to_peptides_psm, select = variables, data.table = FALSE), y = NULL)
      data_list$y = as.numeric(fread(path_to_peptides_psm, select = "PSM Count (unambiguous, <0.01 q-value)", data.table = FALSE)[, 1])
    }
    
    PEPTIDE_DF = do.call("build_peptide_df", data_list)
    rm(data_list)
    protein_df_args = list(protein_name = get_prot_from_EC(PEPTIDE_DF$EC))
  } else if (input_type == "openMS") {
    PEPTIDE_DF = get_peptides_from_idXML(path_to_peptides_psm, PEP)
    PROTEIN_DF_openMS = get_proteins_from_idXML(path_to_peptides_psm)
    protein_name_openMS = get_prot_from_EC(PEPTIDE_DF$EC)
    protein_df_args = list(
      protein_name = PROTEIN_DF_openMS$isoform[match(protein_name_openMS, PROTEIN_DF_openMS$id)],
      id_openMS = PROTEIN_DF_openMS$id[match(protein_name_openMS, PROTEIN_DF_openMS$id)]
    )
  }
  if (!is.null(tpm_path)) {
    protein_df_args$TPM = load_tpm(protein_df_args$protein_name, tpm_path)
  }

  protein_df_args$protein_length = rep(1, length(protein_df_args$protein_name))

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
