normalize_by_gene = function(results_MCMC){
  
  id_genes = as.numeric(factor(results_MCMC$isoform_results$Gene))
  list_unique_gene = unique(id_genes)
  
  id_genes = lapply(list_unique_gene, function(x){which(x == id_genes)})
  norm_gene = lapply(id_genes, function(x){matrix(apply(as.matrix(results_MCMC$PI[, x]), 1, function(y){y/sum(y)}),
                                                  nrow = nrow(results_MCMC$PI), byrow = T)})
  norm_gene = do.call("cbind", norm_gene)
  chain = norm_gene[, order(unlist(id_genes))]
  
  Pi_norm_TPM = lapply(id_genes, function(x){
    gene_tpms = results_MCMC$isoform_results$TPM[x]
    gene_tpms/sum(gene_tpms)
  })
  Pi_norm_TPM = unlist(Pi_norm_TPM)
  
  results_norm_by_gene = data.frame(Pi_norm = colMeans(chain),
                                    Pi_norm_CI_LB = apply(chain, 2, function(x){quantile(x, 0.025)}),
                                    Pi_norm_CI_UB = apply(chain, 2, function(x){quantile(x, 0.975)})
                                    )
  
  cbind(results_MCMC$isoform_results[, c("Gene", "Isoform")], results_norm_by_gene, Pi_norm_TPM)
}

