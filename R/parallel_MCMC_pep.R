parallel_MCMC_pep = function(pep_df, pept_unique_df, prot_df, protein_length, pp, ll, N){
  groups = get_components(pep_df$EC_numeric)
  
  groups_weight = unlist(lapply(groups, function(x){length(x$proteins)}))
  new_order = sort(groups_weight, index.return = TRUE, decreasing = TRUE)$ix
  distributed_index = c()
  for (inc in seq_len(params$n_cores)) {
    distributed_index = c(distributed_index, seq(0, length(new_order), params$n_cores) + inc)
  }
  new_order = na.omit(new_order[distributed_index])
  groups = groups[new_order]
  
  #assign("N_COMPONENTS", length(groups), envir = globalenv())
  components = list_components_for_MCMC_pep(groups, pep_df, pept_unique_df, prot_df, protein_length, pp, ll, N)
  components = aggregate_components_pep(components, ncores = params$n_cores)
  cluster = makeCluster(params$n_cores, type = "PSOCK")
  registerDoParallel(cl = cluster)
  res = foreach(component = iter(components), .combine = cbind) %dorng% {
    names(component) = formalArgs(MCMC_PEP)
    res = do.call(MCMC_PEP, component)
    rbind(res[[1]], res[[2]])
  }
  stopCluster(cluster)
  
  half = nrow(res) / 2
  res = list(PI = res[1:half, ], Y = res[(half+1):nrow(res), ], groups = groups, pept_prot = PEPT_PROT)
  
  res
}