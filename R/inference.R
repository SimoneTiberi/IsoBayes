#' Run latent variable Bayesian model
#'
#' \code{inference} run latent variable Bayesian model taking as input the data created by \code{\link{load_data}}.
#'
#' @param loaded_data \code{list} of \code{data.frame} returned by \code{\link{load_data}}.
#' @param map_iso_gene a character string indicating the path to a csv file with two fields: isoform and gene name. Required to return protein
#' isoforms relative abundances normalized within each gene and to plot results via \code{\link{plot_relative_abundances}}.
#' @param n_cores the number of cores to use during algorithm execution. Default is 1.
#' @param K the number of MCMC iterations. Default is 2000.
#' @param burn_in the number of initial iterations to discard. Default is 1000.
#' @param thin thinning value to apply to the final MCMC chain. Default is 1.
#'
#' @return A \code{list} of two \code{data.frame}: 'isoform_results' and 'normalized_isoform_results' (relative abundances normalized
#' within each gene, if `map_iso_gene` is provided).
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
#'   path_to_tpm = tpm_path)
#'
#' # Define the path to the map_iso_gene.csv file
#' path_to_map_iso_gene = paste0(data_dir, "/map_iso_gene.csv")
#'
#' # Run the algorithm
#' set.seed(169612)
#' results = inference(data_loaded, map_iso_gene = path_to_map_iso_gene)
#' 
#' # results is a list of 2 data.frames:
#' names(results)
#' 
#' # main results:
#' head(results$isoform_results)
#' 
#' # results normalized within genes
#' # (relative abunances add to 1 within each gene):
#' # useful to study alternative splicing within genes:
#' head(results$normalized_isoform_results)
#' 
#' # For more examples see the vignettes:
#' #browseVignettes("SIMBA")
#'
#' @author Simone Tiberi \email{simone.tiberi@unibo.it} and Jordy Bollon \email{jordy.bollon@iit.it}
#'
#' @seealso \code{\link{load_data}} and \code{\link{plot_relative_abundances}}
#'
#' @export
inference = function(loaded_data,
                     map_iso_gene = "",
                     n_cores = 1,
                     K = 2000,
                     burn_in = 1000,
                     thin = 1) {

  input_check_inference(loaded_data, map_iso_gene, n_cores, K, burn_in, thin)

  if (is.null(loaded_data$PROTEIN_DF$TPM)) {
    message("Transcriptomics data not loaded. Inference will be based only on proteomics data.")
    loaded_data$prior = 0
  } else {
    loaded_data$prior = 0.1
  }

  names(loaded_data) = formalArgs(set_MCMC_args)
  args_MCMC = do.call("set_MCMC_args", loaded_data)
  args_MCMC$params = list(n_cores = n_cores, K = K, burn_in = burn_in, thin = thin, PEP = loaded_data$PEP)
  rm(loaded_data)

  if (args_MCMC$params$PEP) {
    results_MCMC = do.call("run_MCMC_pep", args_MCMC)
  } else {
    results_MCMC = do.call("run_MCMC", args_MCMC)
  }

  if (args_MCMC$params$n_cores > 1) {
    old_order = unlist(lapply(results_MCMC$groups, function(x) {x$proteins})
                       )
    old_order = c(old_order, results_MCMC$one_pept_one_prot)
    old_order = sort(old_order, index.return = TRUE)$ix
    results_MCMC$PI = results_MCMC$PI[, old_order]
    results_MCMC$isoform_results = results_MCMC$isoform_results[old_order, ]
  }

  results_MCMC = get_res_MCMC(results_MCMC, args_MCMC$prot_df$protein_name)

  if (!is.null(args_MCMC$prot_df$TPM)) {
    results_MCMC$isoform_results = stat_from_TPM(results_MCMC$isoform_results, args_MCMC$prot_df$TPM, results_MCMC$PI)
    reorder_col = c("Isoform", "Prob_present", "Abundance", "CI_LB", "CI_UB",
                    "Pi", "Pi_CI_LB", "Pi_CI_UB", "TPM", "Log2_FC", "Prob_prot_inc")
  } else {
    reorder_col = c("Isoform", "Prob_present", "Abundance", "CI_LB", "CI_UB",
                    "Pi", "Pi_CI_LB", "Pi_CI_UB")
  }

  if (file.exists(map_iso_gene)) {
    map_iso_gene_file = fread(map_iso_gene, header = FALSE)
    results_MCMC$isoform_results = merge(results_MCMC$isoform_results, map_iso_gene_file, by.x = "Isoform", by.y = "V1")
    colnames(results_MCMC$isoform_results)[ncol(results_MCMC$isoform_results)] = "Gene"
    res_norm = normalize_by_gene(results_MCMC, tpm = !is.null(args_MCMC$prot_df$TPM))
    reorder_col = c("Gene", reorder_col)
  } else {
    res_norm = NULL
  }

  list(isoform_results = results_MCMC$isoform_results[, reorder_col],
       normalized_isoform_results = res_norm
       )
}