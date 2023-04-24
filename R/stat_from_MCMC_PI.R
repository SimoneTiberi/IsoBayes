stat_from_MCMC_PI = function(chain){
  post_mean_probs_iso = colMeans(chain)
  
  # posterior CI for probs_iso:
  CI_probs_iso = apply(results_MCMC$res[[1]], 2, quantile, probs = c(0.025, 0.975))
  CI_probs_iso = matrix(CI_probs_iso, ncol = 2, byrow = TRUE)
  
  res = data.frame(row.names = NULL,
                   post_mean_probs_iso = post_mean_probs_iso,
                   CI_pi_0.025 = CI_probs_iso[,1],
                   CI_pi_0.975 = CI_probs_iso[,2])
  res
}