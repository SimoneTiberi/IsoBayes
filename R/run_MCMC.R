run_MCMC = function(pept_df, prot_df, protein_length, N, M, prior = 0, lib_size, params) {
  set.seed(169612)
  if(prior == 0){
    # vaguely informative prior for "pi" (protein relative abundances)
    # their posterior = dirichlet(Y + 1)
    # where Y = number of protein INTENSITYs
    pp = rep(1, N)
  }else{
    # we assign a small probability to all isoforms, to avoid a 0-prior:
    epsilon = 10^(-5)
    pp = prot_df$TPM/sum(prot_df$TPM) + epsilon
    pp = pp/sum(pp)
    pp = prior * lib_size * pp
  }
  if (params$parallel) {
    res = parallel_MCMC(pept_df, prot_df, protein_length, pp, N, params)
    # normalize PI
    res$PI = t(apply(res$PI, 1, function(x){x/sum(x)}))
  } else {
    res = MCMC(pept_df$EC_numeric, prot_df$Y_unique, protein_length, pp,
               pept_df$Y, N, M, params$K, params$burn_in, params$thin
    )
    # normalize PI
    res$PI = t(apply(res$PI, 1, function(x){x/sum(x)}))
  }
  res
}