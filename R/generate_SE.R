#' Generate SummarizedExperiment object
#'
#' \code{generate_SE} converts the input files, required to run IsoBayes,
#' into a \code{SummarizedExperiment} object.
#' This object should then be passed to \code{\link{input_data}} function.
#' 
#' @param path_to_peptides_psm a character string indicating the path to one of
#' the following files:
#' 
#' i) the psmtsv file from *MetaMorpheus* tool with PSM counts,
#' 
#' ii) the idXML file from *OpenMS* toolkit, or
#' 
#' iii) a \code{data.frame} or a path to a tsv file, formatted as explained
#' in the "Input user-provided data" Section of the vignettes
#' (only when input_type = "other").
#' @param path_to_peptides_intensities (optional) a character string indicating
#' the path to the psmtsv file from *MetaMorpheus* with intensity values.
#' Required if 'abundance_type' equals to "intensities" and
#' input_type equals to "metamorpheus".
#' @param input_type a character string indicating the tool used to obtain
#' the peptides file: "metamorpheus", "openMS" or "other".
#' @param abundance_type a character string indicating the type of input:
#' "psm" or "intensities".
#' @param PEP logical; if TRUE (default), the algorithm will account for
#' the probability that peptides are erroneously detected.
#' If FALSE, PEP is ignored.
#' We suggest using PEP with a weak FDR threshold of 0.1 (default parameters options).
#' This is because peptides with FDR > 0.1 are usually unreliable,
#' and associated to high error probabilities (e.g., PEP > 0.9).
#' @param FDR_thd a numeric value indicating the False Discovery Rate threshold
#' to be used to discard unreliable peptides.
#'
#' @return A \code{SummarizedExperiment} object.
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
#'
#' # For more examples see the vignettes:
#' # browseVignettes("IsoBayes")
#'
#' @author Jordy Bollon \email{jordy.bollon@iit.it}
#' and Simone Tiberi \email{simone.tiberi@unibo.it}
#'
#' @seealso \code{\link{input_data}}
#'
#' @export
generate_SE = function(path_to_peptides_psm = NULL,
                       path_to_peptides_intensities = NULL,
                       input_type = NULL,
                       abundance_type = NULL,
                       PEP = TRUE,
                       FDR_thd = 0.01){
  if (is.null(path_to_peptides_intensities)) {
    path_to_peptides_intensities = ""
  }
  input_check(path_to_peptides_psm, path_to_peptides_intensities,
              input_type, abundance_type, PEP, FDR_thd)
  
  if (input_type == "metamorpheus"){
    SE = from_MM_to_SE(path_to_peptides_psm, path_to_peptides_intensities,
                       abundance_type, PEP, FDR_thd)
  }
  if (input_type == "openMS"){
    SE = from_OpenMS_to_SE(path_to_peptides_psm, PEP, FDR_thd)
  }
  if (input_type == "other"){
    SE = from_UserData_to_SE(path_to_peptides_psm, input_type, PEP, FDR_thd)
  }
  
  SE
}