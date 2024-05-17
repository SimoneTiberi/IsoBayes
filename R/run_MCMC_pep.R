run_MCMC_pep = function(pept_df, pept_unique_df, prot_df, protein_length, N, M,
                        M_unique, prior, lib_size, params) {
  if (prior == 0) {
    # vaguely informative prior for "pi" (protein relative abundances)
    pp = rep(1 / N, N) * N
  } else {
    # we assign a small probability to all isoforms, to avoid a 0-prior:
    epsilon = 10^(-5)
    pp = prot_df$TPM / sum(prot_df$TPM) + epsilon
    pp = pp / sum(pp)
    pp = prior * lib_size * pp
  }
  if (params$n_cores > 1) {
    res = parallel_MCMC_pep(pept_df, pept_unique_df, prot_df,
                            protein_length, pp, N, params)
  } else {
    res = MCMC_PEP(
      pept_df$EC_numeric, pept_df$Y, pept_df$PEP, M,
      pept_unique_df$EC_numeric, pept_unique_df$Y,
      pept_unique_df$PEP, M_unique, protein_length, pp, N, params$K,
      params$burn_in, params$thin
    )
  }
  # normalize PI
  res$PI = t(apply(res$PI, 1, function(x) {
    x / sum(x)
  }))
  res
}