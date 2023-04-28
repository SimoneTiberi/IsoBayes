# take loaded data (list of Peptide and Protein DFs, created at "load" step).
# run MCMC, return results as data.frame
#' @export
inference = function(loaded_data, prior = 0, length_norm = FALSE, parallel = TRUE, n_cores = 2, K = 10000, burn_in = 1000, thin = 5) {
  loaded_data$prior = prior
  loaded_data$length_norm = length_norm

  names(loaded_data) = formalArgs(set_MCMC_args)
  args_MCMC = do.call("set_MCMC_args", loaded_data)
  rm(loaded_data)
  args_MCMC$params = list(parallel = parallel, n_cores = n_cores, K = K, burn_in = burn_in, thin = thin, PEP = loaded_data$PEP)

  if (args_MCMC$params$PEP) {
    results_MCMC = do.call("run_MCMC_pep", args_MCMC)
  } else {
    results_MCMC = do.call("run_MCMC", args_MCMC)
  }

  if (params$parallel) {
    old_order = unlist(lapply(results_MCMC$groups, function(x) {
      x$proteins
    }))
    if (params$PEP) {
      old_order = c(old_order, results_MCMC$pept_prot)
    } else {
      old_order = c(old_order, results_MCMC$one_pept_one_prot)
    }
    old_order = sort(old_order, index.return = T)$ix
    results_MCMC[[1]] = results_MCMC[[1]][, old_order]
    results_MCMC[[2]] = results_MCMC[[2]][, old_order]
  }

  isoform_results = get_res_MCMC(results_MCMC, args_MCMC$prot_df$protein_name, params)
  isoform_results$Y_unique = args_MCMC$prot_df$Y_unique
  isoform_results = stat_from_TPM(isoform_results, args_MCMC$prot_df$TPM, results_MCMC[[1]])

  isoform_results
}
