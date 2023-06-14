parallel_MCMC_pep = function(pep_df, pept_unique_df, prot_df, protein_length, pp, N, params){
  groups = get_components(pep_df$EC_numeric)
  groups = reorder_groups_by_nProteins(groups, params)
  
  components = list_components_for_MCMC_pep(groups, pep_df, pept_unique_df, prot_df, protein_length, pp, N, params)
  one_pept_one_prot = components$one_pept_one_prot
  components = aggregate_components_pep(components$components, ncores = params$n_cores)
  
  cluster = makeCluster(params$n_cores, type = "PSOCK")
  registerDoParallel(cl = cluster)
  
  res = foreach(component = iter(components)) %dorng% {
    names(component) = formalArgs(MCMC_PEP)
    res = do.call(MCMC_PEP, component)
    isoform_results = stat_from_MCMC_Y(res$Y)
    
    list(isoform_results, res$PI)
  }
  stopCluster(cluster)
  
  res = list(PI = do.call("cbind", lapply(res, function(x){x[[2]]})),
             isoform_results = do.call("rbind", lapply(res, function(x){x[[1]]})),
             groups = groups, one_pept_one_prot = one_pept_one_prot
  )
  # normalize PI
  res$PI = t(apply(res$PI, 1, function(x){x/sum(x)}))
  
  res
}
