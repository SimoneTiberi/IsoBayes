get_res_MCMC = function(chains, protein_name, params){
  isoform_results = stat_from_MCMC_Y(chains[[2]])
  isoform_results = cbind(isoform_results, stat_from_MCMC_PI(chains[[1]]))
  isoform_results$proteins = protein_name

  isoform_results
}
