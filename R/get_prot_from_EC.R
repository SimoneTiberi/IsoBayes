get_prot_from_EC = function(EC){
  EC_vector = strsplit(EC, split = "\\|")
  protein_name = sort(unique(unlist(EC_vector)))
  message(glue("{length(protein_name)} protein isoforms."))
  
  protein_name
}
