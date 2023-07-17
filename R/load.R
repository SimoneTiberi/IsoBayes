#' Load and process input data
#'
#' \code{load_data} reads and processes all the input files required to run our model.
#
#' @param path_to_peptides_psm a character string indicating the path to one of the following files:
#' i) the psmtsv file from *MetaMorpheus* tool with PSM counts,
#' ii) the idXML file from *OpenMS* toolkit, or
#' iii) a tsv file, formatted as explained in the "Input user-provided data" Section of the vignettes.
#' For more details on how to create these files see the vignettes.
#' @param path_to_peptides_intensities (optional) a character string indicating the path to the psmtsv file from *MetaMorpheus* with intensity values.
#' Required if 'abundance_type' equals to "intensities".
#' @param path_to_tpm (optional) a character string indicating the path to a tsv file with mRNA isoform TPMs.
#' The tsv file must have 1 row per isoform, and 2 columns:
#' i) 'isoname': a character string indicating the isoform name;
#' ii) 'tpm': a numeric variable indicating the Transcripts Per Million (TPM) count.
#' Column names must be 'isoname' and 'tpm'.
#' @param input_type a character string indicating the tool used to obtain the peptides file: 
#' "metamorpheus" (default), "openMS" or "other".
#' @param abundance_type a character string indicating the type of input: 
#' "psm" (default) or "intensities" (only when "input_type = metamorpheus").
#' @param PEP logical; 
#' if TRUE (default), the algorithm will account for the probability that peptides are erroneously detected.
#' If FALSE, PEP is ignored.
#' Although FDR_thd and PEP can be jointly used; they are meant to be alternatives.
#' @param FDR_thd a numeric value indicating the False Discovery Rate threshold to be used to discard unreliable peptides.
#' Although FDR_thd and PEP can be jointly used; they are meant to be alternatives.
#'
#' @return A \code{list} of \code{data.frame} objects, with the data needed to run \code{\link{inference}} function.
#'
#' @examples
#' # Load internal data to the package:
#' data_dir = system.file("extdata", package = "IsoBayes")
#'
#' # Define the path to the AllPeptides.psmtsv file returned by *MetaMorpheus* tool
#' path_to_peptides_psm = paste0(data_dir, "/AllPeptides.psmtsv")
#'
#' # Load the data
#' data_loaded = load_data(path_to_peptides_psm = path_to_peptides_psm)
#'
#' # For more examples see the vignettes:
#' #browseVignettes("IsoBayes")
#'
#' @author Jordy Bollon \email{jordy.bollon@iit.it} and Simone Tiberi \email{simone.tiberi@unibo.it}
#'
#' @seealso \code{\link{inference}} and \code{\link{plot_relative_abundances}}.
#'
#' @export
load_data = function(path_to_peptides_psm,
                      path_to_peptides_intensities = NULL,
                      path_to_tpm = NULL,
                      input_type = "metamorpheus",
                      abundance_type = "psm",
                      PEP = TRUE,
                      FDR_thd = NULL) {

  if(is.null(FDR_thd)){
    FDR_thd = 1
  }
  if(is.null(path_to_peptides_intensities)){
    path_to_peptides_intensities = ""
  }
  if(is.null(path_to_tpm)){
    path_to_tpm = ""
  }

  input_check(path_to_peptides_psm, path_to_peptides_intensities, path_to_tpm, input_type, abundance_type, PEP, FDR_thd)

  if (input_type == "metamorpheus") {
    variables = c("Protein Accession", "QValue", "Decoy/Contaminant/Target", "Base Sequence")
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
    message("We found:")
    protein_df_args = list(protein_name = get_prot_from_EC(PEPTIDE_DF$EC))
  } else if (input_type == "openMS") {
    PEPTIDE_DF = get_peptides_from_idXML(path_to_peptides_psm, PEP, FDR_thd)
    PROTEIN_DF_openMS = get_proteins_from_idXML(path_to_peptides_psm)
    message("We found:")
    protein_name_openMS = get_prot_from_EC(PEPTIDE_DF$EC)
    protein_df_args = list(protein_name = PROTEIN_DF_openMS$isoform[match(protein_name_openMS, PROTEIN_DF_openMS$id)],
                           id_openMS = PROTEIN_DF_openMS$id[match(protein_name_openMS, PROTEIN_DF_openMS$id)]
                           )
  } else if (input_type == "other") {
    PEPTIDE_DF = as.data.frame(data.table::fread(path_to_peptides_psm))
    if(!is.null(PEPTIDE_DF$FDR)){
      PEPTIDE_DF = PEPTIDE_DF[PEPTIDE_DF$FDR < FDR_thd, ]
    }
    if(PEP & !("PEP" %in% colnames(PEPTIDE_DF))){
      stop(glue("'PEP'=TRUE, but PEP not present in the data. If PEP not available, please set 'PEP'=FALSE."))
    }
    message("We found:")
    protein_df_args = list(protein_name = get_prot_from_EC(PEPTIDE_DF$EC))
  }
  if (path_to_tpm != "") {
    protein_df_args$TPM = load_tpm(protein_df_args$protein_name, path_to_tpm)
  }

  protein_df_args$protein_length = rep(1, length(protein_df_args$protein_name))

  PROTEIN_DF = do.call("data.frame", protein_df_args)
  rm(protein_df_args)

  if (input_type == "metamorpheus") {
    PEPTIDE_DF = PEPTIDE_DF[PEPTIDE_DF$QValue <= FDR_thd, ]
  }

  PEPTIDE_DF = collapse_pept_w_equal_EC(PEPTIDE_DF, PEP)
  PEPTIDE_DF = PEPTIDE_DF[PEPTIDE_DF$Y > 0, ]

  message("After FDR filtering (if used), we will analyze:")
  protein_name_to_keep = get_prot_from_EC(PEPTIDE_DF$EC)
  if (input_type %in% c("metamorpheus", "other")) {
    PROTEIN_DF = PROTEIN_DF[PROTEIN_DF$protein_name %in% protein_name_to_keep, ]
    UNIQUE_PEPT_ABUNDANCE = unique_protein_abundance(PEPTIDE_DF$Y, PEPTIDE_DF$EC, PROTEIN_DF$protein_name)
  } else if (input_type == "openMS") {
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

  # message(glue("Number of multi-mapping peptides: {nrow(PEPTIDE_DF)}"))
  list(PEPTIDE_DF = PEPTIDE_DF, PEPTIDE_DF_unique = PEPTIDE_DF_unique, PROTEIN_DF = PROTEIN_DF, PEP = PEP)
}
