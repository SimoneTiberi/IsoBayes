#' @export
normalize_by_gene = function(results){
  list_unique_gene = unique(results$isoform_results$gene)
  
  for (x in list_unique_gene) {
    sel = x == results$isoform_results$gene
    if(sum(sel)==1){
      results$results_MCMC$PI[, sel] = 1
    }else{
      rel_abundances = results$results_MCMC$PI[, sel]
      results$results_MCMC$PI[, sel] = t(apply(rel_abundances, 1, function(x){x/sum(x)}))
    }
  }
  results_norm_by_gene = data.frame(post_mean_probs_iso = colMeans(results$results_MCMC$PI),
                                    CI_pi_0.025 = apply(results$results_MCMC$PI, 2, function(x){quantile(x, 0.025)}),
                                    CI_pi_0.975 = apply(results$results_MCMC$PI, 2, function(x){quantile(x, 0.975)})
  )
  cbind(results$isoform_results[, c("gene", "Isoform", "probs_TPM")], results_norm_by_gene)
}