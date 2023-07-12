stat_from_MCMC_PI = function(chain){
  # posterior CI for probs_iso:
  CI = hdi(chain, credMass = 0.95)
  
  res = data.frame(row.names = NULL,
                   Pi = colMeans(chain),
                   Pi_CI_LB = CI[1, ],
                   Pi_CI_UB = CI[2, ])
  res
}