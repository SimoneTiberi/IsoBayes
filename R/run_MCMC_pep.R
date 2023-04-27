run_MCMC_pep = function(pept_df, pept_unique_df, prot_df, protein_length, N, M, M_unique, prior = 0, lib_size){
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
      res = parallel_MCMC_pep(pept_df, pept_unique_df, prot_df, protein_length, pp, N)
    }else{
      res = MCMC_PEP(pept_df$EC_numeric, pept_df$Y, pept_df$PEP, M, pept_unique_df$EC_numeric, pept_unique_df$Y,
                     pept_unique_df$PEP, M_unique, protein_length, pp, N, MCMC_CONFIG$K, MCMC_CONFIG$burn_in,
                     MCMC_CONFIG$thin)
    }
  })
  
  list(res = res, time_info = time_info)
}
