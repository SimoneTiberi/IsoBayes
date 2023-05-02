parallel_MCMC = function(pep_df, prot_df, protein_length, pp, N, params){
  groups = get_components(pep_df$EC_numeric)
  groups = reorder_groups_by_nProteins(groups, params)
  components = list_components_for_MCMC(groups, pep_df$EC_numeric, Y_unique = prot_df$Y_unique,
                                        protein_length, pp, pep_df$Y, N, params)
  
  components = aggregate_components(components, ncores = params$n_cores)
  
  one_pept_one_prot = seq_len(N)[-unlist(lapply(groups, function(x){x$proteins}))]
  one_pept_one_prot_Y = prot_df$Y_unique[one_pept_one_prot]
  
  cluster = makeCluster(params$n_cores, type = "PSOCK")
  registerDoParallel(cl = cluster)
  
  res = foreach(component = iter(components), .combine = cbind) %dorng% {
    names(component) = formalArgs(MCMC)
    res = do.call(MCMC, component)
    rbind(res[[1]], res[[2]])
  }
  stopCluster(cluster)
  
  half = nrow(res) / 2
  res = list(PI = res[1:half, ], Y = res[(half+1):nrow(res), ], groups = groups)
  
  final_K = nrow(res$Y)
  res$Y = cbind(res$Y, matrix(rep(one_pept_one_prot_Y, final_K), final_K, byrow = T))
  res$PI = cbind(res$PI, matrix(rep(rep(1, length(one_pept_one_prot_Y)), final_K), final_K, byrow = T))
  res$one_pept_one_prot = one_pept_one_prot
  
  res
}