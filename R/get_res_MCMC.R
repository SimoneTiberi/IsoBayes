get_res_MCMC = function(results_MCMC, protein_name){
  results_MCMC$isoform_results$Isoform = protein_name
  results_MCMC$isoform_results$Gene = gsub("-.*", "", protein_name)
  results_MCMC$isoform_results = cbind(results_MCMC$isoform_results, stat_from_MCMC_PI(results_MCMC$PI))
  results_MCMC$isoform_results$Gene = gsub("-.*", "", results_MCMC$isoform_results$Isoform)

  results_MCMC
}
