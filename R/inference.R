# take loaded data (list of Peptide and Protein DFs, created at "load" step).
# run MCMC, return results as data.frame
#' @export
inference = function(loaded_data, prior = 0.1, length_norm = FALSE, parallel = TRUE, n_cores = 2, K = 10000, burn_in = 1000, thin = 5) {
  if(is.null(loaded_data$PROTEIN_DF$TPM)){
    print("TPM not loaded. Set prior equal to 0.")
    loaded_data$prior = 0
  }else{
    loaded_data$prior = prior
  }
  loaded_data$length_norm = length_norm

  names(loaded_data) = formalArgs(set_MCMC_args)
  args_MCMC = do.call("set_MCMC_args", loaded_data)
  args_MCMC$params = list(parallel = parallel, n_cores = n_cores, K = K, burn_in = burn_in, thin = thin, PEP = loaded_data$PEP)
  rm(loaded_data)

  if (args_MCMC$params$PEP) {
    results_MCMC = do.call("run_MCMC_pep", args_MCMC)
  } else {
    results_MCMC = do.call("run_MCMC", args_MCMC)
  }

  if (args_MCMC$params$parallel) {
    old_order = unlist(lapply(results_MCMC$groups, function(x) {
      x$proteins
    }))
    old_order = c(old_order, results_MCMC$one_pept_one_prot)
    old_order = sort(old_order, index.return = T)$ix
    results_MCMC$PI = results_MCMC$PI[, old_order]
    results_MCMC$isoform_results = results_MCMC$isoform_results[old_order, ]
  }
  
  results_MCMC = get_res_MCMC(results_MCMC, args_MCMC$prot_df$protein_name)
  results_MCMC$isoform_results$Y_unique = args_MCMC$prot_df$Y_unique
  
  if(!is.null(args_MCMC$prot_df$TPM)){
    results_MCMC$isoform_results = stat_from_TPM(results_MCMC$isoform_results, args_MCMC$prot_df$TPM, results_MCMC$PI)
  }
  
  res_norm = normalize_by_gene(results_MCMC)
  
  list(isoform_results = results_MCMC$isoform_results, normalized_isoform_results = res_norm)
}
