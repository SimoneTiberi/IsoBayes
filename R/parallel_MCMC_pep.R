parallel_MCMC_pep = function(pep_df, pept_unique_df, prot_df, protein_length, pp, N, params){
  groups = get_components(pep_df$EC_numeric)
  groups = reorder_groups_by_nProteins(groups, params)
  
  components = list_components_for_MCMC_pep(groups, pep_df, pept_unique_df, prot_df, protein_length, pp, N, params)
  one_pept_one_prot = components$one_pept_one_prot
  components = aggregate_components_pep(components$components, ncores = params$n_cores)
  cluster = makeCluster(params$n_cores, type = "PSOCK")
  registerDoParallel(cl = cluster)
  res = foreach(component = iter(components), .combine = cbind) %dorng% {
    names(component) = formalArgs(MCMC_PEP)
    res = do.call(MCMC_PEP, component)
    rbind(res[[1]], res[[2]])
  }
  stopCluster(cluster)
  half = nrow(res) / 2
  
  list(PI = res[1:half, ], Y = res[(half+1):nrow(res), ], groups = groups, one_pept_one_prot = one_pept_one_prot)
}
