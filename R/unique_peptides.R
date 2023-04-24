unique_peptides = function(EC_numeric){
  EC_length = sapply(EC_numeric, length)
  print("Summary of number of proteins the peptides map to")
  print(summary(EC_length))
  
  sel_unique = EC_length == 1 # peptides with 1 protein only
  print("Percentage of unique peptides")
  mean(sel_unique)
  
  sel_unique
}