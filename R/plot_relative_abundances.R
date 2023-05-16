#' @export
plot_relative_abundances = function(x, gene_id, plot_CI = TRUE, normalize_gene = TRUE){
  sel = x$isoform_results$gene == gene_id
  if(sum(sel) == 1){
    stop(paste0("Only 1 available isoform for gene ", gene_id, ". Plot not returned since the relative abundances is equal to 1."))
  }else if(normalize_gene){
    rel_abundances = x$normalized_isoform_results[sel, "post_mean_probs_iso"]
    CI = x$normalized_isoform_results[sel, c("CI_pi_0.025", "CI_pi_0.975")]
    prop_transc = x$normalized_isoform_results$probs_TPM[sel] / sum(x$normalized_isoform_results$probs_TPM[sel])
  }else{
    rel_abundances = x$isoform_results[sel, "post_mean_probs_iso"]
    CI = x$isoform_results[sel, c("CI_pi_0.025", "CI_pi_0.975")]
    prop_transc = x$isoform_results$probs_TPM[sel]
  }
  # impose an order to the isoforms (according to the over-all relative abudance):
  ord = order(rel_abundances, decreasing = TRUE)
  
  prop_samp = data.frame(feature_id = factor(x$isoform_results$Isoform[sel], levels = x$isoform_results$Isoform[sel][ord]), 
                         prop = c(rel_abundances, prop_transc), #add transcriptomics probs
                         LB = c(pmax(0, CI[, 1]), rep(NA, nrow(CI))), # LB must be >= 0
                         UB = c(pmin(1,  CI[, 2]),  rep(NA, nrow(CI))), # UB must be <= 0
                         Legend = c(rep("Protein", length(x$isoform_results$Isoform[sel])), rep("Transcript", length(x$isoform_results$Isoform[sel])))
                         )
  # Plot the estimated relative abundance:
  ggp = ggplot() +
    geom_bar(data = prop_samp[prop_samp$Legend == "Protein", ], aes_string(x = "feature_id", y = "prop"),
             stat = "identity", position = position_dodge(width = 0.9), colour = "orange", fill = "goldenrod") +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5), 
          axis.text = element_text(size=16), 
          axis.title = element_text(size=14, face="bold"), 
          plot.title = element_text(size=16)
    ) +
    ggtitle(paste0("Gene: ", gene_id)) + xlab("Isoforms") + ylab("Relative abundance") +
    
    geom_point(data = prop_samp, aes_string(x = "feature_id", y = "prop", shape = "Legend", colour = "Legend"), 
               size = 3, alpha = 0.75) +
    scale_color_manual(breaks = c("Protein", "Transcript"), values=c("gold", "brown1"))
  
  if(plot_CI){
    ggp = ggp + geom_errorbar(data = prop_samp, aes_string(x = "feature_id", ymin = "LB", ymax = "UB"),
                              position = position_dodge(width = 0.9), size = 0.5, width = 0.5, alpha = 0.25)
  }
  ggp
}
