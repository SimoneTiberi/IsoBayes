normalize_by_gene = function(results, chains){
  list_unique_gene = unique(results$gene)
  
  for (x in list_unique_gene) {
    sel = x == results$gene
    if(sum(sel)==1){
      chains$PI[, sel] = 1
    }else{
      rel_abundances = chains$PI[, sel]
      chains$PI[, sel] = t(apply(rel_abundances, 1, function(x){x/sum(x)}))
    }
  }
  results_norm_by_gene = data.frame(post_mean_probs_iso = colMeans(chains$PI),
                                    CI_pi_0.025 = apply(chains$PI, 2, function(x){quantile(x, 0.025)}),
                                    CI_pi_0.975 = apply(chains$PI, 2, function(x){quantile(x, 0.975)})
  )
  cbind(results[, c("gene", "Isoform", "probs_TPM")], results_norm_by_gene)
}