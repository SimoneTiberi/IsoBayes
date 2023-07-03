---
bibliography: References.bib
---

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

## Percolator compatibility
SIMBA is compatible with the output returned by *Percolator* tool [@Percolator] provided within the *OpenMS* Toolkit [@OpenMS].
Here, we provide a brief pipeline where several *OpenMS* applications are chained together to generate an idXML file required to run SIMBA with *Percolator* output. The pipeline starts from peptide identification results stored in mzID files.

First, install *OpenMS* toolkit and *Percolator* tool. For instructions on how to install them on your operating system see [OpenMS Installation](https://openms.readthedocs.io/en/latest/openms-applications-and-tools/installation.html) and [Percolator Installation](https://github.com/percolator/percolator).

Next, declare some useful global variable:
``` shell
path_to_data=/path/to/mzIDfiles
path_out=/path/to/output
NTHREADS=4
ENZYME_indexer="Chymotrypsin"
ENZYME_percolator="chymotrypsin"
DECOY_STRING="mz|DECOY_"
fdr=0.01
```

Here we are showing an example with chymotrypsin enzyme. If the data has been generated with another enzyme, please search the corresponding one
inside the following documentations
``` shell
PeptideIndexer --help
PeptideIndexer --help
```
and reset the global variables `ENZYME_indexer` and `ENZYME_percolator` with the correct enzyme.
The pipeline also assume that in the `/path/to/mzIDfiles` folder there is a fasta file listing target and decoy protein isoforms.
The `DECOY_STRING` allow you to change the string needed to identify a decoy in the fasta file.

Now, the next pipeline will create the final idXML suitable for SIMBA.
``` shell
cd $path_out

# convert mzID file into idXML files
for mz in $path_to_data/*.mzID
do
        IDFileConverter -in $mz -threads $NTHREADS -out $mz.idXML
done

# merge the files
IDMerger -in $path_to_data/*.idXML -threads $NTHREADS -merge_proteins_add_PSMs -out $path_out/merge.idXML
rm $path_to_data/*.idXML

# index the peptide file with the fasta file
PeptideIndexer -in $path_out/merge.idXML -enzyme:name $ENZYME_indexer -threads $NTHREADS -decoy_string_position prefix -decoy_string $DECOY_STRING -fasta $path_to_data/genecodeAndDecoy.fasta -out $path_out/merge_index.idXML
rm $path_out/merge.idXML

# run percolator
PercolatorAdapter -in $path_out/merge_index.idXML -enzyme $ENZYME_percolator -threads $NTHREADS -generic_feature_set -score_type pep -out $path_out/merge_index_percolator_pep.idXML
rm $path_out/merge_index.idXML

# Filter out unreliable peptides (FDR greater than $fdr)
FalseDiscoveryRate -in $path_out/merge_index_percolator_pep.idXML -out $path_out/merge_index_percolator_pep_$fdr.idXML -protein false -threads $NTHREADS -FDR:PSM $fdr -algorithm:add_decoy_peptides -algorithm:add_decoy_proteins
rm $path_out/merge_index_percolator_pep.idXML

# Associate each peptite with Posterior Error Probability score
IDScoreSwitcher -in $path_out/merge_index_percolator_pep_$fdr.idXML -out $path_out/merge_index_percolator_pep_switched_$fdr.idXML -new_score 'Posterior Error Probability_score' -new_score_orientation lower_better -new_score_type pep -threads $NTHREADS
rm $path_out/merge_index_percolator_pep_$fdr.idXML
```

For more details on *OpenMS* tools see its [Documentation](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_documentation.html).

# References