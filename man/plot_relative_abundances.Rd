% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_relative_abundances.R
\name{plot_relative_abundances}
\alias{plot_relative_abundances}
\title{Plot isoform results}
\usage{
plot_relative_abundances(
  res_inference,
  gene_id,
  plot_CI = TRUE,
  normalize_gene = TRUE
)
}
\arguments{
\item{res_inference}{\code{list} of two \code{data.frame} objects
returned by \code{\link{inference}}.}

\item{gene_id}{a character string indicating the gene to be plotted.}

\item{plot_CI}{logical; if TRUE (default),
plot 0.95 level Credibility Intervals for each isoform.}

\item{normalize_gene}{logical; if TRUE (default),
plot isoform relative abundances,
normalized within the specified gene (they add to 1 within a gene).}
}
\value{
A \code{ggplot} object, showing isoform relative abundances for a specific gene.
}
\description{
\code{plot_relative_abundances} plots protein isoforms results,
obtained by \code{\link{inference}},
for a specific gene, together with transcripts abundances if available.
}
\examples{
# see the example of inference function:
help(inference)

}
\seealso{
\code{\link{inference}}
}
\author{
Jordy Bollon \email{jordy.bollon@iit.it}
and Simone Tiberi \email{simone.tiberi@unibo.it}
}
