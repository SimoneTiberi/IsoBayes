#' Load and process input data
#'
#' \code{load_data} reads and processes all the input files required to run the model.
#
#' @param path_to_peptides_psm a character string indicating the path to the psmtsv file from *MetaMorpheus* tool,
#' the idXML file from *OpenMS* toolkit or the tsv file with data coming from any bioinformatics tool. For more details on how to create these files
#' see the vignettes.
#' @param path_to_peptides_intensities a character string indicating the path to the psmtsv file from *MetaMorpheus* with intensity values.
#' Required if 'abundance_type' equal to "intensities".
#' @param path_to_tpm a character string indicating the path to a tsv file with mRNA isoform abundances.
#' The tsv file must have two fields for each isoform: 'isoname', a string for the isoform name, and 'tpm' a numeric variable for the corresponing
#' Transcripts Per Million (tpm) count.
#' @param input_type a character string indicating the tool that outputs the peptides file: "metamorpheus", "openMS" or "other",
#' with default "metamorpheus".
#' @param abundance_type a character string indicating the type of input: "psm" or "intensities", with default "psm".
#' @param PEP logical; if FALSE, the algorithm will not take into account the Peptite Error Probability. Default is TRUE.
#' @param FDR_thd a numeric value indicating the False Discovery Rate threshold to be used to discard unreliable peptides.
#'
#' @return A \code{list} with a peptide \code{data.frame}, a peptide \code{data.frame} with unique peptides (only for PEP=TRUE)
#' and a protein \code{data.frame}.
#'
#' @examples
#' # Load internal data to the package:
#' data_dir = system.file("extdata", package = "SIMBA")
#'
#' # Define the path to the AllPeptides.psmtsv file returned by *MetaMorpheus* tool
#' path_to_peptides_psm = paste0(data_dir, "/AllPeptides.psmtsv")
#'
#' # Load the data
#' data_loaded = load_data(path_to_peptides_psm = path_to_peptides_psm)
#'
#' # For more examples see the vignettes:
#' #browseVignettes("SIMBA")
#'
#' @author Simone Tiberi \email{simone.tiberi@unibo.it} and Jordy Bollon \email{jordy.bollon@iit.it}
#'
#' @seealso \code{\link{inference}} and \code{\link{plot_relative_abundances}}.
#'
#' @export
load_data = function(path_to_peptides_psm,
                      path_to_peptides_intensities = "",
                      path_to_tpm = "",
                      input_type = "metamorpheus",
                      abundance_type = "psm",
                      PEP = TRUE,
                      FDR_thd = NULL
) {
  browser()
  if(is.null(FDR_thd)){
    FDR_thd = 1
  }
  input_check(path_to_peptides_psm, path_to_peptides_intensities, path_to_tpm, input_type, abundance_type, PEP, FDR_thd)

  if (input_type == "openMS") {
    abundance_type = "psm"
    message("input_type = 'openMS'. 'FDR_thd' is ignored and abundance_type = 'psm'.")
  }

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
    message("Found:")
    protein_df_args = list(protein_name = get_prot_from_EC(PEPTIDE_DF$EC))
  } else if (input_type == "openMS") {
    PEPTIDE_DF = get_peptides_from_idXML(path_to_peptides_psm, PEP)
    PROTEIN_DF_openMS = get_proteins_from_idXML(path_to_peptides_psm)
    message("Found:")
    protein_name_openMS = get_prot_from_EC(PEPTIDE_DF$EC)
    protein_df_args = list(protein_name = PROTEIN_DF_openMS$isoform[match(protein_name_openMS, PROTEIN_DF_openMS$id)],
                           id_openMS = PROTEIN_DF_openMS$id[match(protein_name_openMS, PROTEIN_DF_openMS$id)]
                           )
  } else if (input_type == "other") {
    PEPTIDE_DF = as.data.frame(data.table::fread(path_to_peptides_psm))
    if(PEP & !("PEP" %in% colnames(PEPTIDE_DF))){
      stop(glue("'PEP'=TRUE, but PEP not present in the data. If PEP not available, please set 'PEP'=FALSE."))
    }
    message("Found:")
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

  message("We may actually detect:")
  protein_name_to_keep = get_prot_from_EC(PEPTIDE_DF$EC)
  if (input_type %in% c("metamorpheus", "other")) {
    PROTEIN_DF = PROTEIN_DF[PROTEIN_DF$protein_name %in% protein_name_to_keep, ]
    UNIQUE_PEPT_ABUNDANCE = unique_protein_abundance(PEPTIDE_DF$Y, PEPTIDE_DF$EC, PROTEIN_DF$protein_name)
  } else if (input_type == "OpenMS") {
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

  message(glue("Number of multi-mapping peptides: {nrow(PEPTIDE_DF)}"))
  list(PEPTIDE_DF = PEPTIDE_DF, PEPTIDE_DF_unique = PEPTIDE_DF_unique, PROTEIN_DF = PROTEIN_DF, PEP = PEP)
}
