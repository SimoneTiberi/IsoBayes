list_components_for_MCMC = function(components, EC_numeric, Y_unique,
                                    protein_length, pp, Y, N, params) {
    vec = rep(NA, N)
    components = lapply(components, function(x) {
        vec[x$proteins] = seq_len(length(x$proteins))
        list(
            EC_numeric = lapply(EC_numeric[x$peptides], function(x) {
                vec[x]
            }),
            Y_unique = Y_unique[x$proteins],
            protein_length = protein_length[x$proteins],
            pp = pp[x$proteins],
            Y = Y[x$peptides],
            N = length(x$proteins),
            M = length(x$peptides),
            K = params$K,
            burn_in = params$burn_in,
            thin = params$thin,
            trace = params$traceplot
        )
    })
    components
}