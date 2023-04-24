stat_from_TPM = function(isoform_results, TPM){
  isoform_results$TPM = TPM
  isoform_results$probs_TPM = TPM/sum(TPM)
  isoform_results$log2_FC = log2(isoform_results$post_mean_probs_iso/isoform_results$probs_TPM)
  
  isoform_results$prob_iso_greater_tpm = sapply(seq_len(nrow(isoform_results)), function(i){
    mean(results_MCMC$res[[1]][, i] > isoform_results$probs_TPM[i])})
  
  isoform_results
}