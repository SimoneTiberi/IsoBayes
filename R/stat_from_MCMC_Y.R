stat_from_MCMC_Y = function(chains){
  # probability of non-zero (i.e., present):
  p_non_zero = apply(chains > 0, 2, mean)
  
  # posterior mean
  post_mean = colMeans(chains)
  
  # 0.95 CI for protein abundance:
  CI = apply(chains, 2, quantile, probs = c(0.025, 0.975))
  CI = matrix(CI, ncol = 2, byrow = TRUE)
  
  out = data.frame(Prob_present = p_non_zero,
                   Abundance = post_mean,
                   CI_LB = CI[,1], 
                   CI_UB = CI[,2]
                   )
  out
}
