#' Load and process input data
#'
#' \code{input_data} reads and processes a \code{SummarizedExperiment} object collecting
#' input data and metadata required to run IsoBayes model.
#
#' @param SE a \code{SummarizedExperiment} object created by \code{\link{generate_SE}} function.
#' Alternatively, this object can be created by the user, following the structure
#' specified in the "Input user-provided data" Section of the vignettes
#' @param path_to_tpm (optional) a \code{data.frame} object or a character string
#' indicating the path to a tsv file with mRNA isoform TPMs.
#' The tsv file must have 1 row per isoform, and 2 columns:
#' i) 'isoname': a character string indicating the isoform name;
#' ii) 'tpm': a numeric variable indicating
#' the Transcripts Per Million (TPM) count.
#' Column names must be 'isoname' and 'tpm'.
#'
#' @return A \code{list} of \code{data.frame} objects, with the data needed
#' to run \code{\link{inference}} function.
#'
#' @examples
#' # Load internal data to the package:
#' data_dir = system.file("extdata", package = "IsoBayes")
#'
#' # Define the path to the AllPeptides.psmtsv file returned by *MetaMorpheus* tool
#' path_to_peptides_psm = paste0(data_dir, "/AllPeptides.psmtsv")
#' 
#' # Generate a SummarizedExperiment object
#' SE = generate_SE(path_to_peptides_psm = path_to_peptides_psm,
#'                  abundance_type = "psm",
#'                  input_type = "metamorpheus"
#'                  )
#' # Load and process SE object
#' data_loaded = input_data(SE)
#'
#' # For more examples see the vignettes:
#' # browseVignettes("IsoBayes")
#'
#' @author Jordy Bollon \email{jordy.bollon@iit.it}
#' and Simone Tiberi \email{simone.tiberi@unibo.it}
#'
#' @seealso \code{\link{generate_SE}}, \code{\link{inference}}
#'
#' @export
input_data = function(SE,
                      path_to_tpm = NULL) {
  if (!is(SE, "SummarizedExperiment")) {
    stop("SE should be a SummarizedExperiment object.")
  }
  if (is.null(metadata(SE)$input_type)){
    metadata(SE)$input_type = "other"
  }
  PROTEIN_DF = data.frame(protein_name = metadata(SE)$protein_name)
  PEPTIDE_DF = data.frame(colData(SE), Y = t(assay(SE)))
  metadata_SE = metadata(SE)
  rm(SE)
  
  if (!is.null(path_to_tpm)) {
    PROTEIN_DF$TPM = load_tpm(PROTEIN_DF$protein_name, path_to_tpm)
  }
  # PROTEIN_DF$protein_length = rep(1, length(PROTEIN_DF$protein_name))
  PROTEIN_DF$id_openMS = metadata_SE$id_openMS
  
  # CALCULATE N_peptides BEFORE COLLAPSING PEPTIDES
  if(metadata_SE$input_type == "metamorpheus"){
    EC = strsplit(PEPTIDE_DF$EC, split = "\\|")
    EC = lapply(EC, unique)
    protein_name_to_keep = sort(unique(unlist(EC)))
    
    if(ncol(PROTEIN_DF) == 1){
      PROTEIN_DF = data.frame(protein_name = PROTEIN_DF[PROTEIN_DF$protein_name %in% protein_name_to_keep,])
    }else{
      PROTEIN_DF = PROTEIN_DF[PROTEIN_DF$protein_name %in% protein_name_to_keep,]
    }
    
    if(metadata_SE$PEP){
      pep_ALL = rep(PEPTIDE_DF$PEP, sapply(EC, length))
      EC = unlist(EC)
      PROTEIN_DF$protein_length = sapply(PROTEIN_DF$protein_name, function(id){
        pep = pep_ALL[EC == id]
        length(pep) - sum(pep)
      })
      rm(pep_ALL)
    }else{
      EC = unlist(EC)
      PROTEIN_DF$protein_length = sapply(PROTEIN_DF$protein_name, function(id){
        sum(EC == id)
      })
    }
    # if length is 0 or NA, set it to 1.
    sel = is.na(PROTEIN_DF$protein_length) | (PROTEIN_DF$protein_length == 0)
    PROTEIN_DF$protein_length[ sel ] = 1
    rm(sel); rm(EC);
  }else{
    PROTEIN_DF$protein_length = 1
  }
  # maybe move protein_length after protein_name_to_keep -> it has to be done before collapsing ECs though
  
  PEPTIDE_DF = collapse_pept_w_equal_EC(PEPTIDE_DF, metadata_SE$PEP)
  PEPTIDE_DF = PEPTIDE_DF[PEPTIDE_DF$Y > 0,]
  
  message("After filtering petides, we will analyze:")
  protein_name_to_keep = get_prot_from_EC(PEPTIDE_DF$EC)
  if (metadata_SE$input_type %in% c("metamorpheus", "other")) {
    PROTEIN_DF = PROTEIN_DF[PROTEIN_DF$protein_name %in% protein_name_to_keep,]
    UNIQUE_PEPT_ABUNDANCE = unique_protein_abundance(PEPTIDE_DF$Y,
                                                     PEPTIDE_DF$EC,
                                                     PROTEIN_DF$protein_name)
  } else if (metadata_SE$input_type == "openMS") {
    PROTEIN_DF = PROTEIN_DF[PROTEIN_DF$id_openMS %in% protein_name_to_keep,]
    UNIQUE_PEPT_ABUNDANCE = unique_protein_abundance(PEPTIDE_DF$Y,
                                                     PEPTIDE_DF$EC,
                                                     PROTEIN_DF$id_openMS)
  }
  PROTEIN_DF$Y_unique = UNIQUE_PEPT_ABUNDANCE$Y_unique
  
  if (metadata_SE$PEP) {
    PEPTIDE_DF_unique = PEPTIDE_DF[UNIQUE_PEPT_ABUNDANCE$sel_unique,]
    PEPTIDE_DF_unique$EC_numeric = UNIQUE_PEPT_ABUNDANCE$EC_numeric[UNIQUE_PEPT_ABUNDANCE$sel_unique]
    PEPTIDE_DF_unique$EC = NULL
  } else {
    PEPTIDE_DF_unique = NULL
  }
  # keep multi-mapping peptides
  PEPTIDE_DF = PEPTIDE_DF[!UNIQUE_PEPT_ABUNDANCE$sel_unique,]
  PEPTIDE_DF$EC_numeric = UNIQUE_PEPT_ABUNDANCE$EC_numeric[!UNIQUE_PEPT_ABUNDANCE$sel_unique]
  PEPTIDE_DF$EC = NULL
  
  overall_abundance = get_overall_abundance(PEPTIDE_DF, PEPTIDE_DF_unique,
                                            PROTEIN_DF$Y_unique)
  
  if (overall_abundance > 2 * 10 ^ 5) {
    PEPTIDE_DF$Y = PEPTIDE_DF$Y / overall_abundance * 10 ^ 5
    # round abundances to closest integer, BUT we add 0.5 so that
    # very small abundances (between 0 and 0.5) are rounded to 1.
    PEPTIDE_DF$Y = round(PEPTIDE_DF$Y + 0.5)
    if (is.null(PEPTIDE_DF_unique)) {
      PROTEIN_DF$Y_unique = PROTEIN_DF$Y_unique / overall_abundance * 10 ^ 5
      PROTEIN_DF$Y_unique = round(PROTEIN_DF$Y_unique + 0.5)
    } else {
      PEPTIDE_DF_unique$Y = PEPTIDE_DF_unique$Y / overall_abundance * 10 ^ 5
      PEPTIDE_DF_unique$Y = round(PEPTIDE_DF_unique$Y + 0.5)
    }
  }
  list(
    PEPTIDE_DF = PEPTIDE_DF,
    PEPTIDE_DF_unique = PEPTIDE_DF_unique,
    PROTEIN_DF = PROTEIN_DF,
    PEP = metadata_SE$PEP
  )
}
