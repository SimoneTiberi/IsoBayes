input_check_plot = function(x, gene_id, plot_CI, normalize_gene){
  if(is.null(x$normalized_isoform_results)){
    stop("List of genes not provided. Plotting results by gene is not possible. Before executing the 'inference' function,
    provide a csv file with the name of protein isoforms and the gene names. For more details see the vignettes.")
  }
  if(all(names(x) != c("isoform_results", "normalized_isoform_results"))){
    stop("Names of 'x' should be: 'isoform_results', 'normalized_isoform_results'")
  }
  if(!is.character(gene_id)){
    stop("gene_id must be a character string indicating the name of a gene.")
  }
  if (!is.logical(plot_CI)) {
    stop("plot_CI must be a boolean value.")
  }
  if (!is.logical(normalize_gene)) {
    stop("normalize_gene must be a boolean value.")
  }
}