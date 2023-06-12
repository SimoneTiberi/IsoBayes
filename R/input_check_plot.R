input_check_plot = function(x, gene_id, plot_CI, normalize_gene){
  if(names(x) != c("isoform_results", "normalized_isoform_results")){
    stop("Names of 'x' should be: 'isoform_results', 'normalized_isoform_results'")
  }
  if(!is.character(gene_id)){
    stop("Input error: gene_id must be a character string indicating the name of a gene.")
  }
  if (!is.logical(plot_CI)) {
    stop("Input error: plot_CI must be a boolean value.")
  }
  if (!is.logical(normalize_gene)) {
    stop("Input error: normalize_gene must be a boolean value.")
  }
}