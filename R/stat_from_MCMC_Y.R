stat_from_MCMC_Y = function(chains){
  # probability of non-zero (i.e., present):
  p_non_zero = apply(chains > 0, 2, mean)
  
  # posterior mean, median and mode for Y abundance
  post_mean = colMeans(chains)
  post_median = apply(chains, 2, quantile, probs = c(0.5))
  post_mode = apply(chains, 2, DescTools::Mode)
  
  # 0.95 CI for protein abundance:
  CI = apply(chains, 2, quantile, probs = c(0.025, 0.975))
  CI = matrix(CI, ncol = 2, byrow = TRUE)
  
  out = data.frame(row.names = NULL,
                   Probability_present = p_non_zero,
                   Posterior_mode = sapply(post_mode, function(y){y[[1]]}),
                   Posterior_median = post_median,
                   Posterior_mean = post_mean,
                   CI_0.025 = CI[,1], 
                   CI_0.975 = CI[,2])
  
  out
}