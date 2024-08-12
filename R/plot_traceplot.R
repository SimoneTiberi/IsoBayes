#' Traceplot of the (thinned) posterior chain of the relative abundance of each protein isoform (i.e., pi).
#'
#' \code{plot_traceplot} plots the traceplot of the (thinned) posterior chain 
#' of the relative abundance of each protein isoform (i.e., pi).
#' The vertical grey dashed line indicates the burn-in 
#' (the iterations on the left side of the burn-in are discarded in posterior analyses).
#' 
#' @param results a \code{list} of \code{\linkS4class{data.frame}} objects, 
#' computed via \code{\link{inference}}.
#' @param protein_id a character, indicating the protein isoform to plot.
#' 
#' @return A \code{gtable} object.
#' 
#' @examples
#' # see the example of inference function:
#' help(inference)
#' 
#' @author Simone Tiberi \email{simone.tiberi@unibo.it}
#' 
#' @seealso \code{\link{inference}}
#' 
#' @export
plot_traceplot = function(results,
                          protein_id){
  if( !is.list(results) ){
    message("'results' must be a 'list' created via 'inference' function.")
    return(NULL)
  }
  if(! "MCMC" %in% names(results)){
    message("'names(results)' must include 'MCMC'.")
    message("Set 'traceplot = TRUE' when running 'inference' function.")
    return(NULL)
  }
  if( is.null(results$MCMC) ){
    message("'results$MCMC' is NULL.")
    message("Set 'traceplot = TRUE' when running 'inference' function.")
    return(NULL)
  }
  
  if(!is.character(protein_id)){
    protein_id = as.character(protein_id)
  }
  if( length(protein_id) != 1){
    message("'protein_id' contains ", length(protein_id), " values: one protein isoform only can be specified.")
    return(NULL)
  }
  
  sel_protein = which( results$MCMC$Isoform == protein_id)
  
  if( length(sel_protein) == 0 ){
    message("'protein_id' not found in 'results$MCMC$Isoform'")
    return(NULL)
  }
  
  n_iter = nrow(results$MCMC$PI)
  burn_in = results$MCMC$thinned_burn_in
  
  DF = data.frame(PI = results$MCMC$PI[,sel_protein],
                  MCMC_iterations = 1:n_iter)
  
  # Plot the estimated average proportions of each groups:
  ggplot() +
    geom_line(data = DF, aes_string(x = "MCMC_iterations", y = "PI")) +
    theme_bw() + 
    theme(axis.text = element_text(size=16), 
          axis.title = element_text(size=14, face="bold"), 
          plot.title = element_text(size=16)) +
    ggtitle(paste("Isoform:", protein_id)) +
    ylab(expression(pi)) +
    geom_vline(xintercept = burn_in, linetype="dashed", 
               colour = "darkgrey")
}