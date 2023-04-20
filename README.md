# SIMBA: Single Isoform protein inference Method via Bayesian Analyses

<img src="inst/extdata/SIMBA.png" width="200" align="right"/> 

`SIMBA` is a Bayesian method to perform inference on single protein isoforms.
Our approach infers the presence/absence of protein isoforms, and also estimates their abundance;
additionally, it provides a measure of the uncertainty of these estimates, via:
i) the posterior probability the a protein isoform is present in the sample;
ii) a posterior credible interval of its abundance.
SIMBA inputs liquid cromatography mass spectrometry (MS) data,
and can work with both PSM counts, and intensities.
When available, trascript isoform abundances (i.e., TPMs) are also incorporated:
TPMs are used to formulate an informative prior for the respective protein isoform relative abundance.
We further identify isoforms where the relative abundance of proteins and transcripts significantly differ.
We use a two-layer latent variable approach to model two sources of uncertainty typical of MS data:
i) peptides may be erroneously detected (even when absent);
ii) many peptides are compatible with multiple protein isoforms.
In the first layer, we sample the presence/absence of each peptide based on its estimated probability 
of being mistakenly detected, also known as PEP.
In the second layer, for peptides that were estimated as being present, 
we allocate their abundance across the protein isoforms they map to.
These two steps allow us to recover the presence and abundance of each protein isoform.

## Bioconductor installation 
`SIMBA` is available on [Bioconductor](https://bioconductor.org/packages/SIMBA) and can be installed with the command:
``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
install.packages("BiocManager")
BiocManager::install("SIMBA")
```

## Vignette
The vignette illustrating how to use the package can be accessed on 
[Bioconductor](https://bioconductor.org/packages/SIMBA)
or from R via:
``` r
vignette("SIMBA")
```
or
``` r
browseVignettes("SIMBA")
```
