#normalize_by_gene = function(results, chain){
#  list_unique_gene = unique(results$gene)
  
#  for (x in list_unique_gene) {
#    sel = x == results$gene
#    if(sum(sel)==1){
#      chain[, sel] = 1
#    }else{
#      rel_abundances = chain[, sel]
#      chain[, sel] = t(apply(rel_abundances, 1, function(x){x/sum(x)}))
#    }
#  }
#  results_norm_by_gene = data.frame(post_mean_probs_iso = colMeans(chain),
#                                    CI_pi_0.025 = apply(chain, 2, function(x){quantile(x, 0.025)}),
#                                    CI_pi_0.975 = apply(chain, 2, function(x){quantile(x, 0.975)})
#  )
#  cbind(results[, c("gene", "Isoform", "probs_TPM")], results_norm_by_gene)
#}

normalize_by_gene = function(results_MCMC){
  id_genes = as.numeric(factor(results_MCMC$isoform_results$gene))
  list_unique_gene = unique(id_genes)
  
  id_genes = lapply(list_unique_gene, function(x){which(x == id_genes)})
  norm_gene = lapply(id_genes, function(x){matrix(apply(as.matrix(results_MCMC$PI[, x]), 1, function(y){y/sum(y)}),
                                                  nrow = nrow(results_MCMC$PI), byrow = T)})
  norm_gene = do.call("cbind", norm_gene)
  
  chain = norm_gene[, order(unlist(id_genes))]
  
  results_norm_by_gene = data.frame(post_mean_probs_iso = colMeans(chain),
                                    CI_pi_0.025 = apply(chain, 2, function(x){quantile(x, 0.025)}),
                                    CI_pi_0.975 = apply(chain, 2, function(x){quantile(x, 0.975)})
                                    )
  
  cbind(results_MCMC$isoform_results[, c("gene", "Isoform", "probs_TPM")], results_norm_by_gene)
}

