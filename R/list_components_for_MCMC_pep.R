list_components_for_MCMC_pep = function(groups, pep_df, pept_unique_df,
                                        prot_df, protein_length, pp, N,
                                        params) {
    vec = rep(NA, N)
    components = lapply(groups, function(x) {
        vec[x$proteins] = seq_len(length(x$proteins))
        pept_unique_in_prot = pept_unique_df$EC_numeric %in% x$proteins

        if (any(pept_unique_in_prot)) {
            EC_numeric_unique = pept_unique_df$EC_numeric[pept_unique_in_prot]
            Y_unique = pept_unique_df$Y[pept_unique_in_prot]
            PEP_unique = pept_unique_df$PEP[pept_unique_in_prot]
        } else {
            EC_numeric_unique = 0
            Y_unique = 0
            PEP_unique = 1
        }

        list(
            EC_numeric = lapply(pep_df$EC_numeric[x$peptides], function(x) {
                vec[x]
            }),
            Y = pep_df$Y[x$peptides],
            PEP = pep_df$PEP[x$peptides],
            M = length(x$peptides),
            EC_numeric_unique = ifelse(EC_numeric_unique == 0, 1,
                                       vec[EC_numeric_unique]),
            Y_unique = Y_unique,
            PEP_unique = PEP_unique,
            M_unique = length(PEP_unique),
            protein_length = protein_length[x$proteins],
            pp = pp[x$proteins],
            N = length(x$proteins),
            K = params$K,
            burn_in = params$burn_in,
            thin = params$thin,
            trace = params$traceplot
        )
    })
    list_pept_prot = get_list_pept_prot(pept_unique_df, groups, protein_length,
                                        pp, params)

    list(components = append(components, list_pept_prot[[1]]),
         one_pept_one_prot = list_pept_prot[[2]])
}