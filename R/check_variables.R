check_variables = function(x){
  if (!("Y" %in% colnames(x))) {
    stop("Peptide abundance 'Y' not found. Check if it is present in the data
          and/or the corresponding column name is 'Y'.")
  }
  if (!is.numeric(x$Y)){
    stop("Peptide abundance 'Y' should be numeric.")
  }
  
  if (!("EC" %in% colnames(x))) {
    stop("Equivalent Classes 'EC' not found. Check if it is present in the data
          and/or the corresponding column name is 'EC'.")
  }
  if (!is.character(x$EC)){
    stop("Equivalent Classes 'EC' should be character.")
  }
  
  if (!is.null(x$FDR)){
    if (!is.numeric(x$FDR)){
      stop("False Discovery Rate 'FDR' should be numeric.")
    }
  }
  
  if (!is.null(x$FDR)){
    if (!is.numeric(x$PEP)){
      stop("Peptide Error Probability 'PEP' should be numeric.")
    }
    if (!("sequence" %in% colnames(x))){
      stop("Peptide name/id/amino acids sequence not found (required
           with PEP). Check if it is present in the data and/or
           the corresponding column name is 'sequence'.")
    }
    if (!is.character(x$sequence)){
      stop("Peptide name/id/amino acids 'sequence' should be character.")
    }
  }
}