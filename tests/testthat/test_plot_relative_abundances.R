test_that("plot_relative_abundances() works faultlessly.", {
    # Load internal data to the package:
    data_dir <- system.file("extdata", package = "IsoBayes")

    # Define the path to the AllPeptides.psmtsv file returned by MetaMorpheus tool
    path_to_peptides_psm <- paste0(data_dir, "/AllPeptides.psmtsv")
    
    # Generate a SummerizedExperiment object
    SE <- generate_SE(path_to_peptides_psm = path_to_peptides_psm,
                      abundance_type = "psm",
                      input_type = "metamorpheus",
                      PEP = FALSE
    )
    # Define the path to the jurkat_isoform_kallisto.tsv with mRNA relative abundance
    tpm_path <- paste0(data_dir, "/jurkat_isoform_kallisto.tsv")

    # Load the data
    data_loaded <- input_data(SE, path_to_tpm = tpm_path)

    # Define the path to the map_iso_gene.csv file
    path_to_map_iso_gene <- paste0(data_dir, "/map_iso_gene.csv")

    # Run the algorithm
    set.seed(169612)
    results <- inference(data_loaded, map_iso_gene = path_to_map_iso_gene,
                         n_cores = 2)

    # Plotting results, NOT normalized
    # (relative abundances add to 1 across all isoforms in the dataset):
    gg <- plot_relative_abundances(results,
        gene_id = "TUBB",
        normalize_gene = FALSE
    )

    expect_is(gg, "ggplot")
})
