#' Run our two-layer latent variable Bayesian model
#'
#' \code{inference} runs our two-layer latent variable Bayesian model,
#' taking as input the data created by \code{\link{input_data}}.
#'
#' @param loaded_data \code{list} of \code{data.frame} objects,
#' returned by \code{\link{input_data}}.
#' @param map_iso_gene (optional) a character string (indicating the path to a csv file with 2 columns),
#'or a data.frame with 2 columns.
#' In both cases, the 1st column must contain the isoform name/id, while the 2nd column has the gene name/id.
#' This argument is required to return protein isoform relative abundances,
#' normalized within each gene
#' (i.e., adding to 1 within a gene), to plot results via \code{\link{plot_relative_abundances}},
#' and to return protein abundances aggregated by gene with HPD credible interval.
#' @param n_cores the number of cores to use during algorithm execution.
#' We suggest increasing the number of threads for large datasets only.
#' @param K the number of MCMC iterations. Minimum 2000.
#' @param burn_in the number of initial iterations to discard. Minimum 1000.
#' @param thin thinning value to apply to the final MCMC chain.
#' Useful for decreasing the memory (RAM) usage.
#' @param traceplot a logical value indicating whether to return the posterior chain
#' of the relative abundances of each protein isoform (i.e., "PI").
#' If TRUE, the posterior chains are stored in 'MCMC' object,
#' and can be plotted via 'plot_traceplot' function.
#'
#' @return A \code{list} of three \code{data.frame} objects: 'isoform_results',
#' and (only if `map_iso_gene` is provided) 'normalized_isoform_results'
#' (relative abundances normalized within each gene)
#' and 'gene_abundance'. For more information about the results stored
#' in the three \code{data.frame} objects, see the vignettes:
#' #browseVignettes("IsoBayes")
#'
#' @examples
#' # Load internal data to the package:
#' data_dir = system.file("extdata", package = "IsoBayes")
#'
#' # Define the path to the AllPeptides.psmtsv file returned by MetaMorpheus tool
#' path_to_peptides_psm = paste0(data_dir, "/AllPeptides.psmtsv")
#' 
#' # Generate a SummarizedExperiment object
#' SE = generate_SE(path_to_peptides_psm = path_to_peptides_psm,
#'                  abundance_type = "psm",
#'                  input_type = "metamorpheus"
#'                  )
#' # Define the path to the jurkat_isoform_kallisto.tsv with mRNA relative abundance
#' tpm_path = paste0(data_dir, "/jurkat_isoform_kallisto.tsv")
#'            
#' # Load and process SE object
#' data_loaded = input_data(SE, path_to_tpm = tpm_path)
#'
#' # Define the path to the map_iso_gene.csv file.
#' # Alternatively a data.frame can be used (see documentation).
#' path_to_map_iso_gene = paste0(data_dir, "/map_iso_gene.csv")
#'
#' # Run the algorithm
#' set.seed(169612)
#' results = inference(data_loaded, map_iso_gene = path_to_map_iso_gene, traceplot = TRUE)
#'
#' # Results is a list of 3 data.frames:
#' names(results)
#'
#' # Main results:
#' head(results$isoform_results)
#'
#' # Results normalized within genes
#' # (relative abunances add to 1 within each gene):
#' # useful to study alternative splicing within genes:
#' head(results$normalized_isoform_results)
#'
#' # Gene abundance
#' head(results$gene_abundance)
#'
#' # results normalized within genes (total abundance of each gene),
#' # useful to study alternative splicing within genes:
#' head(results$normalized_isoform_results)
#'
#' # Plotting results, normalizing within genes
#' # (relative abundances add to 1 within each gene):
#' plot_relative_abundances(results,
#'     gene_id = "TUBB",
#'     normalize_gene = TRUE)
#'
#' # Plotting results, NOT normalized
#' # (relative abundances add to 1 across all isoforms in the dataset):
#' plot_relative_abundances(results,
#'     gene_id = "TUBB",
#'     normalize_gene = FALSE)
#' 
#' # Visualize MCMC chain for isoforms "TUBB-205", "TUBB-206", and  "TUBB-208"
#' # To visualize traceplots, set "traceplot" to TRUE when running "inference" function
#' plot_traceplot(results, "TUBB-205")
#' plot_traceplot(results, "TUBB-206")
#' plot_traceplot(results, "TUBB-208")
#'
#' # For more examples see the vignettes:
#' # browseVignettes("IsoBayes")
#'
#' @author Jordy Bollon \email{jordy.bollon@iit.it}
#' and Simone Tiberi \email{simone.tiberi@unibo.it}
#'
#' @seealso \code{\link{input_data}} and \code{\link{plot_relative_abundances}}
#'
#' @export
inference = function(loaded_data,
                     map_iso_gene = NULL,
                     n_cores = 1,
                     K = 2000,
                     burn_in = 1000,
                     thin = 1,
                     traceplot = FALSE) {
  if(!is.logical(traceplot)){
    message("'traceplot' must be 'TRUE' or 'FALSE'.")
    return(NULL)
  }
  if( !(traceplot %in% c(TRUE,FALSE)) ){
    message("'traceplot' must be 'TRUE' or 'FALSE'.")
    return(NULL)
  }
  
  if (is.null(map_iso_gene)) {
    map_iso_gene = ""
  }
  
  input_check_inference(loaded_data, map_iso_gene, n_cores, K, burn_in, thin)
  
  if (is.null(loaded_data$PROTEIN_DF$TPM)) {
    message("Transcriptomics data not loaded.
                Inference will be based only on proteomics data.")
    loaded_data$prior = 0
  } else {
    loaded_data$prior = 0.1
  }
  
  map_iso_gene_file = NULL
  if ( is.data.frame(map_iso_gene) ) {
    if(ncol(map_iso_gene) != 2){
      stop("'map_iso_gene' must be a character string to a csv file with 2 columns, 
           or directly a data.frame with 2 columns: in both cases, with isoform (1st column) and gene (2nd column) ids.")
    }
    map_iso_gene_file = map_iso_gene
    colnames(map_iso_gene_file) = c("V1", "V2")
    rm(map_iso_gene)
  }else{
    if ( is.character(map_iso_gene) ){
      map_iso_gene_file = fread(map_iso_gene, header = FALSE)
      if(ncol(map_iso_gene_file) != 2){
        stop("'map_iso_gene' must be a character string to a csv file with 2 columns, 
           or directly a data.frame with 2 columns: in both cases, with isoform (1st column) and gene (2nd column) ids.")
      }
    }
  }

  names(loaded_data) = formalArgs(set_MCMC_args)
  args_MCMC = do.call("set_MCMC_args", loaded_data)
  args_MCMC$params = list(n_cores = n_cores, K = K, burn_in = burn_in,
                          thin = thin, PEP = loaded_data$PEP,
                          traceplot = traceplot)
  sel_unique = loaded_data$prot_df$Y_unique > 0
  rm(loaded_data)
  
  if (args_MCMC$params$PEP){
    results_MCMC = do.call("run_MCMC_pep", args_MCMC)
  }else{
    results_MCMC = do.call("run_MCMC", args_MCMC)
  }
  
  if (args_MCMC$params$n_cores > 1){
    old_order = unlist(lapply(results_MCMC$groups, function(x){
      x$proteins
    }))
    old_order = c(old_order, results_MCMC$one_pept_one_prot)
    old_order = sort(old_order, index.return = TRUE)$ix
    results_MCMC$PI = results_MCMC$PI[, old_order]
    results_MCMC$Y = results_MCMC$Y[, old_order]
    if(traceplot){
      results_MCMC$PI_burn_in = results_MCMC$PI_burn_in[, old_order]
    }
  }
  
  if(traceplot){
    MCMC = list( PI = rbind(results_MCMC$PI_burn_in, results_MCMC$PI),
                 Isoform = args_MCMC$prot_df$protein_name,
                 thinned_burn_in = burn_in/thin)
  }else{
    MCMC = NULL
  }
  results_MCMC$PI_burn_in = NULL
  
  results_MCMC$isoform_results = stat_from_MCMC_Y(results_MCMC$Y)
  results_MCMC = get_res_MCMC(results_MCMC, args_MCMC$prot_df$protein_name)
  
  if (!is.null(args_MCMC$prot_df$TPM)) {
    results_MCMC$isoform_results = stat_from_TPM(results_MCMC$isoform_results,
                                                 args_MCMC$prot_df$TPM,
                                                 results_MCMC$PI)
    reorder_col = c(
      "Isoform", "Prob_present", "Abundance", "CI_LB", "CI_UB",
      "Pi", "Pi_CI_LB", "Pi_CI_UB", "TPM", "Log2_FC", "Prob_prot_inc"
    )
  } else {
    reorder_col = c(
      "Isoform", "Prob_present", "Abundance", "CI_LB", "CI_UB",
      "Pi", "Pi_CI_LB", "Pi_CI_UB"
    )
  }
  
  if ( !is.null(map_iso_gene_file) ) {
    results_MCMC = map_isoform_to_gene(results_MCMC, map_iso_gene_file)
    
    res_norm = normalize_by_gene(results_MCMC, tpm = !is.null(args_MCMC$prot_df$TPM))
    reorder_col = c("Gene", reorder_col)
    
    results_MCMC$Y = aggregate.data.frame(t(results_MCMC$Y),
                                          by = list(results_MCMC$isoform_results$Gene),
                                          FUN = sum)
    Gene = results_MCMC$Y[, 1]
    results_MCMC$Y = t(results_MCMC$Y[, -1])
    # 0.95 CI for protein abundance:
    CI = hdi(results_MCMC$Y, credMass = 0.95)
    gene_abundance = data.frame(
      Gene = Gene,
      Abundance = colMeans(results_MCMC$Y),
      CI_LB = CI[1, ],
      CI_UB = CI[2, ]
    )
  } else {
    res_norm = NULL
    gene_abundance = NULL
  }
  if (!args_MCMC$params$PEP) {
    # add a small threshold to isoform with unique peptides
    results_MCMC$isoform_results$Prob_present[sel_unique] = results_MCMC$isoform_results$Prob_present[sel_unique] + 0.01
  }
  
  list(
    isoform_results = results_MCMC$isoform_results[, reorder_col],
    normalized_isoform_results = res_norm,
    gene_abundance = gene_abundance,
    MCMC = MCMC
  )
}
