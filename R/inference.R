# take loaded data (list of Peptide and Protein DFs, created at "load" step).
# run MCMC, return results as data.frame
#' @export
inference = function(loaded_data, prior = 0, parallel = TRUE, K = 10000, burn_in = 1000, thin = 5){
  MCMC_CONFIG = list(parallel = parallel, K = K, burn_in = burn_in, thin = thin, PEP = loaded_data$PEP)
  
  loaded_data$prior = prior
  names(loaded_data) = formalArgs(set_MCMC_args)
  args_MCMC = do.call("set_MCMC_args", loaded_data)
  rm(loaded_data)
  
  if(MCMC_CONFIG$PEP){
    results_MCMC = do.call("run_MCMC_pep", args_MCMC)
  }else{
    results_MCMC = do.call("run_MCMC", args_MCMC)
  }
  
  isoform_results = get_res_MCMC(results_MCMC$res, args_MCMC$prot_df$protein_name)
  isoform_results$Y_unique = args_MCMC$prot_df$Y_unique
  isoform_results = stat_from_TPM(isoform_results, args_MCMC$prot_df$TPM)
  
  isoform_results
}