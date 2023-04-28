list_components_for_MCMC = function(components, EC_numeric, Y_unique, protein_length, pp, Y, ll, N){
  vec = rep(NA, N)
  components = lapply(components, function(x){
    vec[x$proteins] = seq_len(length(x$proteins))
    list(EC_numeric = lapply(EC_numeric[x$peptides], function(x){vec[x]}),
         Y_unique = Y_unique[x$proteins],
         protein_length = protein_length[x$proteins],
         pp = pp[x$proteins],
         Y = Y[x$peptides],
         N = length(x$proteins),
         M = length(x$peptides),
         K = CONFIG$MCMC$K,
         burn_in = CONFIG$MCMC$burn_in,
         thin = CONFIG$MCMC$thin,
         ll = ll)
  })
  components
}