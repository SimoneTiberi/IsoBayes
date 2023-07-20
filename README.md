# IsoBayes: Isoform-level Bayesian proteogenomics inference

<img src="inst/extdata/IsoBayes.png" width="200" align="right"/> 

`IsoBayes` is a Bayesian method to perform inference on single protein isoforms.
Our approach infers the presence/absence of protein isoforms, and also estimates their abundance;
additionally, it provides a measure of the uncertainty of these estimates, via:
i) the posterior probability that a protein isoform is present in the sample;
ii) a posterior credible interval of its abundance.
`IsoBayes` inputs liquid cromatography mass spectrometry (MS) data,
and can work with both PSM counts, and intensities.
When available, trascript isoform abundances (i.e., TPMs) are also incorporated:
TPMs are used to formulate an informative prior for the respective protein isoform relative abundance.
We further identify isoforms where the relative abundance of proteins and transcripts significantly differ.
We use a two-layer latent variable approach to model two sources of uncertainty typical of MS data:
i) peptides may be erroneously detected (even when absent);
ii) many peptides are compatible with multiple protein isoforms.
In the first layer, we sample the presence/absence of each peptide based on its estimated probability 
of being mistakenly detected, also known as PEP (i.e., posterior error probability).
In the second layer, for peptides that were estimated as being present, 
we allocate their abundance across the protein isoforms they map to.
These two steps allow us to recover the presence and abundance of each protein isoform.

## Installation 
`IsoBayes` should appear on Bioconductor in a few weeks.
In the meantime, you can install it from GitHub:
``` r
devtools::install_github("SimoneTiberi/IsoBayes")
```


## Vignette
The vignette illustrating how to use the package can be obtained by running the `IsoBayes.Rmd` file in the vignettes folder.

## Input data
*IsoBayes* works with the output of both *MetaMorpheus* (MM), or *Percolator* (via the *OpenMS* toolkit).
Other tools can be used as well to generate out input data, as long as the input files follow the structure mentioned in the "Input user-provided data" Section of the vignettes.
In our benchmarks, we tested our model using both *MetaMorpheus* and *Percolator* data, and obtained slightly better results, and a shorter runtime with *MetaMorpheus*.

### *MetaMorpheus* pipeline
To generate the MM output files required to run *IsoBayes*, we need to execute the following commands:

* Install *MetaMorpheus* via [*Conda*](https://docs.conda.io/en/latest/miniconda.html):
```shell
conda install -c conda-forge metamorpheus
```

* Inside the folder with the configuration (.toml), spectra (.mzML or .raw) and database (.xml) files run:
```shell
metamorpheus -t Task1-SearchTaskconfig.toml Task2-CalibrateTaskconfig.toml Task3-SearchTaskconfig.toml Task4-GPTMDTaskconfig.toml Task5-SearchTaskconfig.toml -s 04-30-13_CAST_Frac4_6uL.raw 04-30-13_CAST_Frac5_4uL.raw -d uniprot-mouse-reviewed-1-24-2018.xml.gz uniprot-cRAP-1-24-2018.xml.gz
```
or
```shell
metamorpheus -t Task1-SearchTaskconfig.toml Task2-CalibrateTaskconfig.toml Task3-SearchTaskconfig.toml Task4-GPTMDTaskconfig.toml Task5-SearchTaskconfig.toml -s mzML/04-30-13_CAST_Frac4_6uL.mzML mzML/04-30-13_CAST_Frac5_4uL.mzML -d uniprot-mouse-reviewed-1-24-2018.xml.gz uniprot-cRAP-1-24-2018.xml.gz
```
There are several ways to install and run MM. For more details see the MM [tutorial](https://github.com/smith-chem-wisc/MetaMorpheus/wiki/Getting-Started#test-installation-via-net-core-dll---linux-macos-windows), where you can also find the example files used here.

### *Percolator* pipeline
We provide a brief pipeline where several *OpenMS* applications are chained together to generate an idXML file required to run IsoBayes with *Percolator* output. The pipeline starts from peptide identification results stored in mzID files.

First, install *OpenMS* toolkit and *Percolator* tool.
For instructions on how to install them on your operating system see [OpenMS Installation](https://openms.readthedocs.io/en/latest/openms-applications-and-tools/installation.html) and [Percolator Installation](https://github.com/percolator/percolator).

Next, declare some useful global variable:
``` shell
path_to_data=/path/to/mzIDfiles
path_out=/path/to/output
NTHREADS=4
ENZYME_indexer="Chymotrypsin"
ENZYME_percolator="chymotrypsin"
DECOY_STRING="mz|DECOY_"
fdr=1
```

Below, we show an example with chymotrypsin enzyme.
If the data was generated with another enzyme, please search for the corresponding enzyme in the following documentation below, and reset the global variables `ENZYME_indexer` and `ENZYME_percolator` with the correct enzyme.
``` shell
PeptideIndexer --help
PercolatorAdapter --help
```

This pipeline also assumes that in the `/path/to/mzIDfiles` folder there is a fasta file listing target and decoy protein isoforms.
The `DECOY_STRING` allows you to change the string needed to identify a decoy in the fasta file.

``` shell
cd $path_out

# convert mzID files into idXML files
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

# Estimate the false discovery rate on peptide level using decoy searches and keep the ones with FDR < $fdr
FalseDiscoveryRate -in $path_out/merge_index_percolator_pep.idXML -out $path_out/merge_index_percolator_pep_$fdr.idXML -protein false -threads $NTHREADS -FDR:PSM $fdr -algorithm:add_decoy_peptides -algorithm:add_decoy_proteins
rm $path_out/merge_index_percolator_pep.idXML

# Associate each peptite with Posterior Error Probability score
IDScoreSwitcher -in $path_out/merge_index_percolator_pep_$fdr.idXML -out $path_out/merge_index_percolator_pep_switched_$fdr.idXML -new_score 'Posterior Error Probability_score' -new_score_orientation lower_better -new_score_type pep -threads $NTHREADS
rm $path_out/merge_index_percolator_pep_$fdr.idXML
```

For more details on *OpenMS* tools see its [Documentation](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/nightly/html/TOPP_documentation.html).

# References

Röst, H. L., Sachsenberg, T., Aiche, S., Bielow, C., Weisser, H., Aicheler, F., ... & Kohlbacher, O. (2016). OpenMS: a flexible open-source software platform for mass spectrometry data analysis. *Nature methods*, 13(9), 741-748.

The, M., MacCoss, M. J., Noble, W. S., & Käll, L. (2016). Fast and accurate protein false discovery rates on large-scale proteomics data sets with percolator 3.0. *Journal of the American Society for Mass Spectrometry*, 27, 1719-1727.

Solntsev, S. K., Shortreed, M. R., Frey, B. L., & Smith, L. M. (2018). Enhanced global post-translational modification discovery with MetaMorpheus. *Journal of proteome research*, 17(5), 1844-1851.
