get_res_MCMC = function(chains, protein_name){
  if(MCMC_CONFIG$parallel){
    vec = 1:ncol(chains[[1]])
    old_order = unlist(lapply(chains$groups, function(x){x$proteins}))
    if(MCMC_CONFIG$PEP){
      old_order = c(old_order, chains$pept_prot)
    }else{
      old_order = c(old_order, chains$one_pept_one_prot)
    }
    old_order = sort(old_order, index.return = T)$ix
    chains[[1]] = chains[[1]][, old_order]
    chains[[2]] = chains[[2]][, old_order]
  }
  isoform_results = stat_from_MCMC_Y(chains[[2]])
  isoform_results = cbind(isoform_results, stat_from_MCMC_PI(chains[[1]]))
  isoform_results$proteins = protein_name

  isoform_results
}
