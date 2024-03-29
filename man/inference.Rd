% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inference.R
\name{inference}
\alias{inference}
\title{Run our two-layer latent variable Bayesian model}
\usage{
inference(
  loaded_data,
  map_iso_gene = NULL,
  n_cores = 1,
  K = 2000,
  burn_in = 1000,
  thin = 1
)
}
\arguments{
\item{loaded_data}{\code{list} of \code{data.frame} objects,
returned by \code{\link{input_data}}.}

\item{map_iso_gene}{(optional) a character string indicating
the path to a csv file with two columns:
the 1st one containing the isoform id, and the 2nd one with the gene name.
This argument is required to return protein isoform relative abundances,
normalized within each gene
(i.e., adding to 1 within a gene), to plot results via \code{\link{plot_relative_abundances}},
and to return protein abundances aggregated
by gene with HPD credible interval.}

\item{n_cores}{the number of cores to use during algorithm execution.
We suggest increasing the number of threads for large datasets only.}

\item{K}{the number of MCMC iterations. Minimum 2000.}

\item{burn_in}{the number of initial iterations to discard. Minimum 1000.}

\item{thin}{thinning value to apply to the final MCMC chain.
Useful for decreasing the memory (RAM) usage.}
}
\value{
A \code{list} of three \code{data.frame} objects: 'isoform_results',
and (only if `map_iso_gene` is provided) 'normalized_isoform_results'
(relative abundances normalized within each gene)
and 'gene_abundance'. For more information about the results stored
in the three \code{data.frame} objects, see the vignettes:
#browseVignettes("IsoBayes")
}
\description{
\code{inference} runs our two-layer latent variable Bayesian model,
taking as input the data created by \code{\link{input_data}}.
}
\examples{
# Load internal data to the package:
data_dir = system.file("extdata", package = "IsoBayes")

# Define the path to the AllPeptides.psmtsv file returned by MetaMorpheus tool
path_to_peptides_psm = paste0(data_dir, "/AllPeptides.psmtsv")

# Generate a SummarizedExperiment object
SE = generate_SE(path_to_peptides_psm = path_to_peptides_psm,
                 abundance_type = "psm",
                 input_type = "metamorpheus"
                 )
# Define the path to the jurkat_isoform_kallisto.tsv with mRNA relative abundance
tpm_path = paste0(data_dir, "/jurkat_isoform_kallisto.tsv")
           
# Load and process SE object
data_loaded = input_data(SE, path_to_tpm = tpm_path)

# Define the path to the map_iso_gene.csv file
path_to_map_iso_gene = paste0(data_dir, "/map_iso_gene.csv")

# Run the algorithm
set.seed(169612)
results = inference(data_loaded, map_iso_gene = path_to_map_iso_gene)

# Results is a list of 3 data.frames:
names(results)

# Main results:
head(results$isoform_results)

# Results normalized within genes
# (relative abunances add to 1 within each gene):
# useful to study alternative splicing within genes:
head(results$normalized_isoform_results)

# Gene abundance
head(results$gene_abundance)

# For more examples see the vignettes:
# browseVignettes("IsoBayes")

}
\seealso{
\code{\link{input_data}} and \code{\link{plot_relative_abundances}}
}
\author{
Jordy Bollon \email{jordy.bollon@iit.it}
and Simone Tiberi \email{simone.tiberi@unibo.it}
}
