get_res_MCMC = function(chains, protein_name){
  isoform_results = stat_from_MCMC_Y(chains[[2]], protein_name)
  isoform_results = cbind(isoform_results, stat_from_MCMC_PI(chains[[1]]))
  isoform_results$gene = gsub("-.*", "", protein_name)

  isoform_results
}
