test_that("load_data() works faultlessly.", {
  # Load internal data to the package:
  data_dir = system.file("extdata", package = "SIMBA")
  
  # Define the path to the AllPeptides.psmtsv file returned by *MetaMorpheus* tool
  path_to_peptides_psm = paste0(data_dir, "/AllPeptides.psmtsv")
  
  # Load the data
  data_loaded = load_data(path_to_peptides_psm = path_to_peptides_psm)
  
  expect_is(data_loaded, "list")
  expect_true(length(data_loaded) == 4)
})