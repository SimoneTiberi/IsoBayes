check_and_load_tpm = function(tpm) {
  if (is.character(tpm)) {
    if (!file.exists(tpm)) {
      stop("Path to tpm file does not exist!")
    } else {
      tpm = fread(file = tpm, header = TRUE)
    }
  } else if (!is.data.frame(tpm)) {
    stop(
      "'tpm' should be a data.frame object or a character path to a tsv file with isoforms name and tpm values."
    )
  }
  if (!all(c("isoname", "tpm") %in% colnames(tpm))) {
    stop("Columns names of the tpm file should be: 'isoname' and 'tpm'.")
  }
  
  tpm
}