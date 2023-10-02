test_that("inference() works faultlessly.", {
    # Load internal data to the package:
    data_dir <- system.file("extdata", package = "IsoBayes")

    # Define the path to the AllPeptides.psmtsv file returned by MetaMorpheus tool
    path_to_peptides_psm <- paste0(data_dir, "/AllPeptides.psmtsv")

    # Define the path to the jurkat_isoform_kallisto.tsv with mRNA relative abundance
    tpm_path <- paste0(data_dir, "/jurkat_isoform_kallisto.tsv")

    # Load the data
    data_loaded <- load_data(
        path_to_peptides_psm = path_to_peptides_psm,
        path_to_tpm = tpm_path
    )

    # Define the path to the map_iso_gene.csv file
    path_to_map_iso_gene <- paste0(data_dir, "/map_iso_gene.csv")

    # Run the algorithm
    set.seed(169612)
    results <- inference(data_loaded, map_iso_gene = path_to_map_iso_gene)

    expect_is(results, "list")
    expect_true(length(results) == 3)
    expect_is(results[[1]], "data.frame")
})
