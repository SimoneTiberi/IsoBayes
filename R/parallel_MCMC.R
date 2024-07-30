parallel_MCMC = function(pep_df, prot_df, protein_length, pp, N, params) {
    groups = get_components(pep_df$EC_numeric)
    groups = reorder_groups_by_nProteins(groups, params)
    components = list_components_for_MCMC(groups, pep_df$EC_numeric,
        Y_unique = prot_df$Y_unique,
        protein_length, pp, pep_df$Y, N, params
    )

    components = aggregate_components(components, ncores = params$n_cores)

    one_pept_one_prot = seq_len(N)[-unlist(lapply(groups, function(x) {
        x$proteins
    }))]
    one_pept_one_prot_Y = prot_df$Y_unique[one_pept_one_prot]

    cluster = makeCluster(params$n_cores, type = "PSOCK")
    registerDoParallel(cl = cluster)

    res = foreach(component = iter(components)) %dorng% {
        names(component) = formalArgs(MCMC)
        res = do.call(MCMC, component)
        res
    }
    stopCluster(cluster)

    PI_one_pept_one_prot = MCMC_Unique(
        one_pept_one_prot_Y, pp[one_pept_one_prot],
        length(one_pept_one_prot_Y), params$K,
        params$burn_in, params$thin, params$traceplot
    )
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
        groups = groups
    )
    res$PI = cbind(res$PI, PI_one_pept_one_prot$PI)
    res$Y = cbind(res$Y, matrix(rep(one_pept_one_prot_Y,
                                    round((params$K - params$burn_in) / params$thin)),
        nrow(res$Y),
        byrow = TRUE
    ))
    res$PI_burn_in = cbind(res$PI_burn_in, PI_one_pept_one_prot$PI_burn_in)
    res$one_pept_one_prot = one_pept_one_prot

    res
}