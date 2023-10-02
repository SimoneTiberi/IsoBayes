stat_from_TPM = function(isoform_results, TPM, chain) {
    isoform_results$TPM = TPM
    P_TPM = TPM / sum(TPM)
    isoform_results$Log2_FC = log2(isoform_results$Pi / P_TPM)

    isoform_results$Prob_prot_inc = vapply(seq_len(nrow(isoform_results)),
                                           function(i) {
        mean(chain[, i] > P_TPM[i])
    }, FUN.VALUE = numeric(1))

    isoform_results
}