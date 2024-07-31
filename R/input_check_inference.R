input_check_inference = function(loaded_data, map_iso_gene, n_cores, K,
                                 burn_in, thin) {
  if (all(names(loaded_data) != c(
    "PEPTIDE_DF", "PEPTIDE_DF_unique",
    "PROTEIN_DF", "PEP"
  ))) {
    stop("Names of 'loaded_data' should be: 'PEPTIDE_DF',
        'PEPTIDE_DF_unique', 'PROTEIN_DF', 'PEP'")
  }
  if ( (!is.character(map_iso_gene)) & (!is.data.frame(map_iso_gene)) ) {
    stop("'map_iso_gene' must be a character string, or a data.frame.")
  }
  if ( is.character(map_iso_gene) ) {
    if (map_iso_gene == "") {
      message("'map_iso_gene' not specified.
            We return results without gene normalization.")
    } else if (!file.exists(map_iso_gene)) {
      message("'map_iso_gene' does not exist.
            We return results without gene normalization.")
    }
  }
  if (n_cores != round(n_cores) || n_cores < 0) {
    stop("'n_cores' must be an integer > 0.")
  }
  if (burn_in != round(burn_in) || burn_in < 1000) {
    stop("'burn_in' should be an integer >= 1000.")
  }
  if (K != round(K) || K < 2000) {
    stop("'K' should be an integer >= 2000.")
  }
  if (thin != round(thin) || thin < 1) {
    stop("'thin' should be an integer >= 1.")
  }
  if (round((K - burn_in) / thin) < 1000) {
    stop(glue("With K={K} iterations, burn_in={burn_in} and thin={thin}
    the final chains will have less than 1000 points.
              Increase K or decrease the burn_in or the thin value."))
  }
}