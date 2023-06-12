input_check_inference = function(loaded_data, prior, parallel, n_cores, K, burn_in, thin){
  if(all(names(loaded_data) != c("PEPTIDE_DF", "PEPTIDE_DF_unique", "PROTEIN_DF", "PEP"))){
    stop("Names of 'loaded_data' should be: 'PEPTIDE_DF', 'PEPTIDE_DF_unique', 'PROTEIN_DF', 'PEP'")
  }
  if (prior < 0 & prior > 1) {
    stop("Input error: prior must be a numeric value between 0 and 1.")
  }
  if (!is.logical(parallel)) {
    stop("Input error: parallel must be a boolean value.")
  }
  if(n_cores != round(n_cores) || n_cores < 0){
    stop("Input error: n_cores must be an integer > 0.")
  }
  if(burn_in != round(burn_in) || burn_in < 1000){
    stop("Input error: burn_in should be an integer >= 1000.")
  }
  if(K != round(K) || K < 2000){
    stop("Input error: K should be an integer >= 2000.")
  }
  if(thin != round(thin) || thin < 1){
    stop("Input error: thin should be an integer >= 1.")
  }
  if(round((K-burn_in)/thin) < 1000){
    stop(glue("With K={K} iterations, burn_in={burn_in} and thin={thin} the final chains will have less than 1000 points.
                Increase K or decrease the burn_in or the thin value."))
  }
}