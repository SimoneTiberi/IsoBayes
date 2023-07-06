#' Plot isoform results
#'
#' \code{plot_relative_abundances} plots protein isoforms results for a specific gene together with transcripts abundances if available.
#'
#' @param res_inference \code{list} of \code{data.frame} returned by \code{\link{inference}}.
#' @param gene_id a character string indicating the gene to be plot.
#' @param plot_CI logical; if TRUE (default), plot Credibility Intervals for each isoform.
#' @param normalize_gene logical; if TRUE (default), plot isoforms relative abundances normalized within the specified gene.
#'
#' @return A plot showing isoforms relative abundances for a specific gene.
#'
#' @examples
#' # Load internal data to the package:
#' data_dir = system.file("extdata", package = "SIMBA")
#'
#' # Define the path to the AllPeptides.psmtsv file returned by MetaMorpheus tool
#' path_to_peptides_psm = paste0(data_dir, "/AllPeptides.psmtsv")
#'
#' # Define the path to the jurkat_isoform_kallisto.tsv with mRNA relative abundance
#' tpm_path = paste0(data_dir, "/jurkat_isoform_kallisto.tsv")
#'
#' # Load the data
#' data_loaded = load_data(
#'   path_to_peptides_psm = path_to_peptides_psm,
#'   path_to_tpm = tpm_path
#' )
#'
#' # Define the path to the map_iso_gene.csv file
#' path_to_map_iso_gene = paste0(data_dir, "/map_iso_gene.csv")
#'
#' # Run the algorithm
#' set.seed(169612)
#' results = inference(data_loaded, map_iso_gene = path_to_map_iso_gene)
#'
#' # Plotting results
#' plot_relative_abundances(results, gene_id = "HLA")
#' 
#' # For more examples see the vignettes:
#' #browseVignettes("SIMBA")
#'
#' @author Simone Tiberi \email{simone.tiberi@unibo.it} and Jordy Bollon \email{jordy.bollon@iit.it}
#'
#' @seealso \code{\link{load_data}}, \code{\link{inference}}
#'
#' @export
plot_relative_abundances = function(res_inference,
                                    gene_id,
                                    plot_CI = TRUE,
                                    normalize_gene = TRUE) {
  
  input_check_plot(res_inference, gene_id, plot_CI, normalize_gene)

  sel = res_inference$isoform_results$Gene == gene_id
  if (sum(sel) == 1) {
    stop(paste0("Only 1 available isoform for gene ", gene_id, ". Plot not returned since the relative abundances is equal to 1."))
  } else if (normalize_gene) {
    df_sub = res_inference$normalized_isoform_results[sel, ]
    rel_abundances = df_sub$Pi_norm
    CI = df_sub[, c("Pi_norm_CI_LB", "Pi_norm_CI_UB")]
    prop_transc = df_sub$Pi_norm_TPM
  } else {
    df_sub = res_inference$isoform_results[sel, ]
    rel_abundances = df_sub$Pi
    CI = df_sub[, c("Pi_CI_LB", "Pi_CI_UB")]
    prop_transc = df_sub$TPM / sum(res_inference$isoform_results$TPM)
  }
  # impose an order to the isoforms (according to the overall relative abudances)
  ord = order(rel_abundances, decreasing = TRUE)

  prop_samp = data.frame(feature_id = factor(res_inference$isoform_results$Isoform[sel],
                                             levels = res_inference$isoform_results$Isoform[sel][ord]
                                             ),
                         prop = c(rel_abundances, prop_transc), # add transcriptomics probs
                         LB = c(pmax(0, CI[, 1]), rep(NA, length(prop_transc))), # LB must be >= 0
                         UB = c(pmin(1, CI[, 2]), rep(NA, length(prop_transc))), # UB must be <= 0
                         Legend = c(rep("Protein", length(rel_abundances)), rep("Transcript", length(prop_transc)))
                         )
  # Plot the estimated relative abundances:
  ggp = ggplot() +
    geom_bar(data = prop_samp[prop_samp$Legend == "Protein", ], aes_string(x = "feature_id", y = "prop"),
             stat = "identity", position = position_dodge(width = 0.9), colour = "orange", fill = "goldenrod"
             ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_text(size = 16)
          ) +
    ggtitle(paste0("Gene: ", gene_id)) +
    xlab("Isoforms") +
    ylab("Relative abundance") +
    geom_point(data = prop_samp, aes_string(x = "feature_id", y = "prop", shape = "Legend", colour = "Legend"),
               size = 3, alpha = 0.75
               ) +
    scale_color_manual(breaks = c("Protein", "Transcript"), values = c("gold", "brown1"))

  if (plot_CI) {
    ggp = ggp + geom_errorbar(data = prop_samp, aes_string(x = "feature_id", ymin = "LB", ymax = "UB"),
                              position = position_dodge(width = 0.9), size = 0.5, width = 0.5, alpha = 0.25
                              )
  }
  ggp
}
