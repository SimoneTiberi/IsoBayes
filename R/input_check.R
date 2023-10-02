input_check = function(path_to_peptides_psm, path_to_peptides_intensities,
                        path_to_tpm, input_type, abundance_type, PEP, FDR_thd) {
    character_inputs = c(
        path_to_peptides_intensities, path_to_tpm, input_type,
        abundance_type
    )

    args_name = c(
        "path_to_peptides_intensities", "path_to_tpm", "input_type",
        "abundance_type"
    )

    for (i in 1:length(character_inputs)) {
        if (!is(character_inputs[i], "character")) {
            stop(glue("{args_name[i]} must be a character string."))
        }
    }

    if (!(input_type %in% c("openMS", "metamorpheus", "other"))) {
        stop("Invalid input_type. Choose one of 'openMS',
             'metamorpheus' or 'other'.")
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

    if (input_type == "other") {
        if (!is.character(path_to_peptides_psm) & !is.data.frame(path_to_peptides_psm)) {
            stop("path_to_peptides_psm must be a character string or
                 a data.frame.")
        }
    } else {
        if (!is.character(path_to_peptides_psm)) {
            stop("path_to_peptides_psm must be a character string.")
        }
    }

    if (input_type == "openMS" & abundance_type == "intensities") {
        stop("With input_type = 'openMS' only psm abundance can be inputted.
         Please set abundance_type equal to 'psm'.")
    }

    if (!is.data.frame(path_to_peptides_psm)) {
        vec_files = c(
            path_to_peptides_psm, path_to_peptides_intensities,
            path_to_tpm
        )
        args_name = c(
            "path_to_peptides_psm", "path_to_peptides_intensities",
            "path_to_tpm"
        )
    } else {
        vec_files = c(path_to_peptides_intensities, path_to_tpm)
        args_name = c("path_to_peptides_intensities", "path_to_tpm")
    }
    file_exist = file.exists(vec_files)
    file_not_specified = c(vec_files) == ""
    check_path = (file_exist + file_not_specified) == 1
    if (!all(check_path)) {
        first_error_arg = args_name[which(!check_path)[1]]
        stop(glue("{first_error_arg} does not exist."))
    }
}