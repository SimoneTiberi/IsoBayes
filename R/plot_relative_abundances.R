#' Plot isoform results
#'
#' \code{plot_relative_abundances} plot isoform results for a specific gene.
#'
#' @param x list of dataframes returned by \code{\link{inference}}.
#' @param gene_id a character string indicating the gene to be plot.
#' @param plot_CI a boolean value to plot Credibility Interval.
#' @param normalize_gene a boolean value to plot results normalized by gene.
#'
#' @return A plot with relative frequencies of isoforms detected for a specific gene.
#' @examples
#' @author name
#'
#' @seealso links
#'
#' @export
plot_relative_abundances = function(x, gene_id, plot_CI = TRUE, normalize_gene = TRUE){

  input_check_plot(x, gene_id, plot_CI, normalize_gene)
  
  sel = x$isoform_results$Gene == gene_id
  if(sum(sel) == 1){
    stop(paste0("Only 1 available isoform for gene ", gene_id, ". Plot not returned since the relative abundances is equal to 1."))
  }else if(normalize_gene){
    df_sub = x$normalized_isoform_results[sel, ]
    rel_abundances = df_sub$Pi_norm
    CI = df_sub[, c("Pi_norm_CI_LB", "Pi_norm_CI_UB")]
    prop_transc = df_sub$Pi_norm_TPM
  }else{
    df_sub = x$isoform_results[sel, ]
    rel_abundances = df_sub$Pi
    CI = df_sub[, c("Pi_CI_LB", "Pi_CI_UB")]
    prop_transc = df_sub$TPM / sum(x$isoform_results$TPM)
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
