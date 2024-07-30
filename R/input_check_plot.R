input_check_plot = function(x, gene_id, plot_CI, normalize_gene) {
    if (!is(x, "list")) {
        stop("'res_inference' should be a list of data.frame objects created by function 'inference'.")
    }
    if (!all(c("isoform_results", "normalized_isoform_results") %in% names(x))) {
        stop("Names of 'res_inference' should contain: 'isoform_results', and
             'normalized_isoform_results'")
    }
    if (is.null(x$normalized_isoform_results)) {
        stop("Data.frame with genes and isoform names not provided.
        Plotting results by gene is not possible.
        Before executing the 'inference' function, provide a csv file with
        the name of protein isoforms and the gene names.
        For more details see the vignettes.")
    }
    if (!is.character(gene_id)) {
        stop("gene_id must be a character string indicating
             the name of a gene.")
    }
    if (!is.logical(plot_CI)) {
        stop("plot_CI must be a boolean value.")
    }
    if (!is.logical(normalize_gene)) {
        stop("normalize_gene must be a boolean value.")
    }
}