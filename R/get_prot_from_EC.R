get_prot_from_EC = function(EC){
  EC_vector = strsplit(EC, split = "\\|")
  protein_name = sort(unique(unlist(EC_vector)))
  print(glue("total number of proteins: {length(protein_name)}"))
  
  protein_name
}