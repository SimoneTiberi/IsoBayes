unique_peptides = function(EC_numeric){
  EC_length = vapply(EC_numeric, length, FUN.VALUE = numeric(1) )
  #message("Summary of number of proteins the peptides map to:")
  #message(summary(EC_length))
  
  sel_unique = EC_length == 1 # peptides with 1 protein only
  message(glue("Percentage of unique peptides: {round(mean(sel_unique), 4)*100}%"))
  
  sel_unique
}