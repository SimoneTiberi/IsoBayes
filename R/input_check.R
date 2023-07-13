input_check = function(path_to_peptides_psm, path_to_peptides_intensities,
                       path_to_tpm, input_type, abundance_type, PEP, FDR_thd) {
  character_inputs = c(
    class(path_to_peptides_psm), class(path_to_peptides_intensities),
    class(path_to_tpm), class(input_type), class(abundance_type)
  )

  if (!all(character_inputs == "character")) {
    args_name = c("path_to_peptides_psm", "path_to_peptides_intensities", "path_to_tpm", "input_type", "abundance_type")
    first_error_arg = args_name[which(character_inputs != "character")[1]]
    stop(glue("{first_error_arg} must be a character string."))
  }
  
  if (input_type == "openMS" & abundance_type == 'intensities') {
    stop("With input_type = 'openMS' only psm abundance can be inputted. Please set abundance_type equal to 'psm'.")
  }
  
  file_exist = file.exists(path_to_peptides_psm, path_to_peptides_intensities, path_to_tpm)
  file_not_specified = c(path_to_peptides_psm, path_to_peptides_intensities, path_to_tpm) == ""
  check_path = (file_exist + file_not_specified) == 1
  if(!all(check_path) ){
    args_name = c("path_to_peptides_psm", "path_to_peptides_intensities", "path_to_tpm")
    first_error_arg = args_name[which(!check_path)[1]]
    stop(glue("{first_error_arg} does not exist."))
  }
  
  if (!(input_type %in% c("openMS", "metamorpheus", "other"))) {
    stop("Invalid input_type. Choose one of 'openMS', 'metamorpheus' or 'other'.")
  }
  if (!(abundance_type %in% c("psm", "intensities"))) {
    stop("Invalid input_type. Choose one of 'psm' or 'intensities'.")
  }
  if (!is.logical(PEP)) {
    stop("PEP must be a boolean value.")
  }
  if (FDR_thd < 0 | FDR_thd > 1) {
    stop("FDR_thd must be a numeric value between 0 and 1.")
    }
}
