input_check = function(path_to_peptides_psm, path_to_peptides_intensities, tpm_path, input_type, abundance_type, PEP, FDR_thd) {
  character_inputs = c(
    class(path_to_peptides_psm), class(path_to_peptides_intensities),
    class(tpm_path), class(input_type), class(abundance_type)
  )

  if (!all(character_inputs == "character")) {
    args_name = c("path_to_peptides_psm", "path_to_peptides_intensities", "tpm_path", "input_type", "abundance_type")
    first_error_arg = args_name[which(character_inputs != "character")[1]]
    stop(glue("Input error: {first_error_arg} must be a character string."))
  }
  if (!(input_type %in% c("openMS", "metamorpheus"))) {
    stop("Invalid input_type. Choose one of 'openMS' or 'metamorpheus'.")
  }
  if (!(abundance_type %in% c("psm", "intensities"))) {
    stop("Invalid input_type. Choose one of 'psm' or 'intensities'.")
  }
  if (!is.logical(PEP)) {
    stop("Input error: PEP must be a boolean value.")
  }
  if (FDR_thd < 0 & FDR_thd > 1) {
    stop("Input error: FDR_thd must be a numeric value between 0 and 1.")
  }
}