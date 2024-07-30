parallel_MCMC_pep = function(pep_df, pept_unique_df, prot_df, protein_length,
                             pp, N, params) {
    groups = get_components(pep_df$EC_numeric)
    groups = reorder_groups_by_nProteins(groups, params)

    components = list_components_for_MCMC_pep(groups, pep_df, pept_unique_df,
                                              prot_df, protein_length,
                                              pp, N, params)
    one_pept_one_prot = components$one_pept_one_prot
    components = aggregate_components_pep(components$components,
                                          ncores = params$n_cores)

    cluster = makeCluster(params$n_cores, type = "PSOCK")
    registerDoParallel(cl = cluster)

    res = foreach(component = iter(components)) %dorng% {
        names(component) = formalArgs(MCMC_PEP)
        res = do.call(MCMC_PEP, component)
        res
    }
    stopCluster(cluster)

    res = list(
        PI = do.call("cbind", lapply(res, function(x) {
            x$PI
        })),
        Y = do.call("cbind", lapply(res, function(x) {
            x$Y
        })),
        PI_burn_in = do.call("cbind", lapply(res, function(x) {
          x$PI_burn_in
        })),
        groups = groups, one_pept_one_prot = one_pept_one_prot
    )
    res
}