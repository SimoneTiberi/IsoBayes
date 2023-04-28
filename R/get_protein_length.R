get_protein_length = function(protein_name, path_fasta){
  protein_length = fasta.seqlengths(path_fasta)
  names(protein_length) = gsub("mz\\|", "", names(protein_length))
  names(protein_length) = gsub("\\|.*", "", names(protein_length))
  
  protein_length = protein_length[names(protein_length) %in% protein_name]
  
  protein_length[match(protein_name, names(protein_length))]
}