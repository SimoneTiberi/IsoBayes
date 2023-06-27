#' Run latent variable Bayesian model
#'
#' \code{inference} run latent variable Bayesian model
#'
#' @param loaded_data list of \code data.frame} returned by \code{\link{load}}.
#' @param prior a numeric value indicating how much weight to assign on trascriptomic prior information.
#' @param parallel a boolean value to enable parallel computation.
#' @param n_cores number of cores to be used during parallel computation.
#' @param K number of MCMC iterations. # modificabile dall'utente?
#' @param burn_in number of initial iterations to discard. # modificabile dall'utente?
#' @param thin thinning value to apply to the final MCMC chain. # modificabile dall'utente?
#'
#' @return A list of two dataframe: 'isoform_results' and "normalized_isoform_results" (results normalized by gene).
#' @examples
#' @author name
#'
#' @seealso links
#'
#' @export
inference = function(loaded_data, prior = 0.1, map_iso_gene = "",  parallel = FALSE, n_cores = 2, K = 10000, burn_in = 1000, thin = 5) {

  input_check_inference(loaded_data, prior, map_iso_gene, parallel, n_cores, K, burn_in, thin)
  
  if(is.null(loaded_data$PROTEIN_DF$TPM)){
    print("TPM not loaded. Set prior equal to 0.")
    loaded_data$prior = 0
  }else{
    loaded_data$prior = prior
  }

  names(loaded_data) = formalArgs(set_MCMC_args)
  args_MCMC = do.call("set_MCMC_args", loaded_data)
  args_MCMC$params = list(parallel = parallel, n_cores = n_cores, K = K, burn_in = burn_in, thin = thin, PEP = loaded_data$PEP)
  rm(loaded_data)

  if (args_MCMC$params$PEP) {
    results_MCMC = do.call("run_MCMC_pep", args_MCMC)
  } else {
    results_MCMC = do.call("run_MCMC", args_MCMC)
  }

  if (args_MCMC$params$parallel) {
    old_order = unlist(lapply(results_MCMC$groups, function(x) {
      x$proteins
    }))
    old_order = c(old_order, results_MCMC$one_pept_one_prot)
    old_order = sort(old_order, index.return = T)$ix
    results_MCMC$PI = results_MCMC$PI[, old_order]
    results_MCMC$isoform_results = results_MCMC$isoform_results[old_order, ]
  }
  
  results_MCMC = get_res_MCMC(results_MCMC, args_MCMC$prot_df$protein_name)
  
  if(!is.null(args_MCMC$prot_df$TPM)){
    results_MCMC$isoform_results = stat_from_TPM(results_MCMC$isoform_results, args_MCMC$prot_df$TPM, results_MCMC$PI)
  }
  
  reorder_col = c("Isoform", "Prob_present", "Abundance", "CI_LB", "CI_UB",
                  "Pi", "Pi_CI_LB", "Pi_CI_UB", "TPM", "Log2_FC", "Prob_prot_inc"
                  )
  
  if(file.exists(map_iso_gene)){
    map_iso_gene_file = fread(map_iso_gene, header = FALSE)
    results_MCMC$isoform_results = merge(results_MCMC$isoform_results, map_iso_gene_file, by.x = "Isoform", by.y = "V1")
    colnames(results_MCMC$isoform_results)[ncol(results_MCMC$isoform_results)] = "Gene"
    res_norm = normalize_by_gene(results_MCMC)
    reorder_col = c("Gene", reorder_col)
  }else{
    res_norm = NULL
  }
  
  list(isoform_results = results_MCMC$isoform_results[, reorder_col],
       normalized_isoform_results = res_norm)
}