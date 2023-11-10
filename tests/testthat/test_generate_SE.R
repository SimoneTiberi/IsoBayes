test_that("generate_SE() works faultlessly.", {
  # Load internal data to the package:
  data_dir <- system.file("extdata", package = "IsoBayes")
  
  # Define the path to the AllPeptides.psmtsv file returned by *MetaMorpheus* tool
  path_to_peptides_psm <- paste0(data_dir, "/AllPeptides.psmtsv")
  
  # Generate a SummerizedExperiment object
  SE = generate_SE(path_to_peptides_psm = path_to_peptides_psm,
                   abundance_type = "psm",
                   input_type = "metamorpheus"
                   )
 
  expect_is(SE, "SummarizedExperiment")
})
