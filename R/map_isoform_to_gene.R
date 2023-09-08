map_isoform_to_gene = function(results_MCMC, map_iso_gene_file){
  
  results_MCMC$isoform_results$id = seq_len(nrow(results_MCMC$isoform_results))
  results_MCMC$isoform_results = merge(results_MCMC$isoform_results, 
                                       map_iso_gene_file, by.x = "Isoform", by.y = "V1")
  colnames(results_MCMC$isoform_results)[ncol(results_MCMC$isoform_results)] = "Gene"
  
  if(ncol(results_MCMC$PI) > nrow(results_MCMC$isoform_results)){
    iso_not_found = ncol(results_MCMC$PI) - nrow(results_MCMC$isoform_results)
    message(glue("In map_iso_gene file, {iso_not_found} isoforms were not found and they will be removed from results normalized within each gene."))
  }
  results_MCMC$PI = results_MCMC$PI[, results_MCMC$isoform_results$id]
  results_MCMC$Y = results_MCMC$Y[, results_MCMC$isoform_results$id]
  
  results_MCMC
}
