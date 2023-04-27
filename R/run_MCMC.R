run_MCMC = function(pept_df, prot_df, protein_length, N, M, prior = 0, lib_size){
  set.seed(169612)
  time_info = system.time({
    if(prior == 0){
      # vaguely informative prior for "pi" (protein relative abundances)
      # their posterior = dirichlet(Y + 1)
      # where Y = number of protein INTENSITYs
      pp = rep(1/N, N)
    }else{
      # we assign a small probability to all isoforms, to avoid a 0-prior:
      epsilon = 10^(-5)
      pp = prot_df$TPM/sum(prot_df$TPM) + epsilon
      pp = pp/sum(pp)
      pp = prior * lib_size * pp
    }
    if(FALSE){
      res = parallel_MCMC(pept_df, prot_df, protein_length, pp, N)
    }else{
      res = MCMC(pept_df$EC_numeric, prot_df$Y_unique, protein_length, pp, pept_df$Y, N, M,
                 MCMC_CONFIG$K, MCMC_CONFIG$burn_in, MCMC_CONFIG$thin)
    }
  })
  
  list(res = res, time_info = time_info)
}
