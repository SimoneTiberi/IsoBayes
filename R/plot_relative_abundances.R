#' Plot isoform results
#'
#' \code{plot_relative_abundances} plots protein isoforms results,
#' obtained by \code{\link{inference}},
#' for a specific gene, together with transcripts abundances if available.
#'
#' @param res_inference \code{list} of two \code{data.frame} objects
#' returned by \code{\link{inference}}.
#' @param gene_id a character string indicating the gene to be plotted.
#' @param plot_CI logical; if TRUE (default),
#' plot 0.95 level Credibility Intervals for each isoform.
#' @param normalize_gene logical; if TRUE (default),
#' plot isoform relative abundances,
#' normalized within the specified gene (they add to 1 within a gene).
#'
#' @return A \code{ggplot} object, showing isoform relative abundances for a specific gene.
#'
#' @examples
#' # see the example of DifferentialRegulation function:
#' help(inference)
#'
#' @author Jordy Bollon \email{jordy.bollon@iit.it}
#' and Simone Tiberi \email{simone.tiberi@unibo.it}
#'
#' @seealso \code{\link{inference}}
#'
#' @export
plot_relative_abundances = function(res_inference,
                                     gene_id,
                                     plot_CI = TRUE,
                                     normalize_gene = TRUE) {
    input_check_plot(res_inference, gene_id, plot_CI, normalize_gene)

    sel = res_inference$isoform_results$Gene == gene_id
    if (sum(sel) == 1) {
        stop("Only 1 available isoform for gene ", gene_id,
             ". Plot not returned since the relative abundances is equal to 1.")
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

    prop_samp = data.frame(
        feature_id = factor(res_inference$isoform_results$Isoform[sel],
            levels = res_inference$isoform_results$Isoform[sel][ord]
        ),
        prop = c(rel_abundances, prop_transc), 
        LB = c(pmax(0, CI[, 1]), rep(NA, length(prop_transc))), 
        UB = c(pmin(1, CI[, 2]), rep(NA, length(prop_transc))), 
        Legend = c(rep("Protein", length(rel_abundances)),
                   rep("Transcript", length(prop_transc)))
    )
    # Plot the estimated relative abundances:
    ggp = ggplot() +
        geom_bar(
            data = prop_samp[prop_samp$Legend == "Protein", ],
            aes_string(x = "feature_id", y = "prop"),
            stat = "identity", position = position_dodge(width = 0.9),
            colour = "orange", fill = "goldenrod"
        ) +
        theme_bw() +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 0.5),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 14, face = "bold"),
            plot.title = element_text(size = 16)
        ) +
        ggtitle(paste0("Gene: ", gene_id)) +
        xlab("Isoforms") +
        ylab("Relative abundance") +
        geom_point(
            data = prop_samp, aes_string(x = "feature_id", y = "prop",
                                         shape = "Legend", colour = "Legend"),
            size = 3, alpha = 0.75
        ) +
        scale_color_manual(breaks = c("Protein", "Transcript"),
                           values = c("gold", "brown1"))

    if (plot_CI) {
        ggp = ggp + geom_errorbar(
            data = prop_samp, aes_string(x = "feature_id",
                                         ymin = "LB", ymax = "UB"),
            position = position_dodge(width = 0.9), size = 0.5,
            width = 0.5, alpha = 0.25
        )
    }
    ggp
}