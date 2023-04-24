convert_EC_to_num = function(EC, vec_names){
  # each element corresponds to the row of proteins/gene in PROTEIN_DF
  EC_numeric = lapply(EC, function(X){match(X, vec_names)})
  EC_numeric = lapply(EC_numeric, unique)  # remove potential duplicated proteins in the EC
  
  # filter eventual NA's (there shouldn't be any!)
  sel_NA = is.na(EC_numeric)
  EC_numeric[!sel_NA]
}