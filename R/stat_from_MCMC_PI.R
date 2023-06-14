stat_from_MCMC_PI = function(chain){
  # posterior CI for probs_iso:
  CI = apply(chain, 2, quantile, probs = c(0.025, 0.975))
  CI = matrix(CI, ncol = 2, byrow = TRUE)
  
  res = data.frame(row.names = NULL,
                   Pi = colMeans(chain),
                   Pi_CI_LB = CI[,1],
                   Pi_CI_UB = CI[,2])
  res
}