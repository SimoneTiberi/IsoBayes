normalize_by_gene = function(results_MCMC, tpm) {
    id_genes = as.numeric(factor(results_MCMC$isoform_results$Gene))
    list_unique_gene = unique(id_genes)

    id_genes = lapply(list_unique_gene, function(x) {
        which(x == id_genes)
    })
    norm_gene = lapply(id_genes, function(x) {
        matrix(
            apply(as.matrix(results_MCMC$PI[, x]), 1, function(y) {
                y / sum(y)
            }),
            nrow = nrow(results_MCMC$PI), byrow = TRUE
        )
    })
    norm_gene = do.call("cbind", norm_gene)
    chain = norm_gene[, order(unlist(id_genes))]

    if (tpm) {
        Pi_norm_TPM = lapply(id_genes, function(x) {
            gene_tpms = results_MCMC$isoform_results$TPM[x]
            gene_tpms / sum(gene_tpms)
        })
        Pi_norm_TPM = unlist(Pi_norm_TPM)
    } else {
        Pi_norm_TPM = NULL
    }

    CI = hdi(chain, credMass = 0.95)
    results_norm_by_gene = data.frame(
        Pi_norm = colMeans(chain),
        Pi_norm_CI_LB = CI[1, ],
        Pi_norm_CI_UB = CI[2, ]
    )
    if (tpm) {
        cbind(results_MCMC$isoform_results[, c("Gene", "Isoform")],
              results_norm_by_gene, Pi_norm_TPM)
    } else {
        cbind(results_MCMC$isoform_results[, c("Gene", "Isoform")],
              results_norm_by_gene)
    }
}