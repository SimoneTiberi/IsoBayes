aggregate_sim_pep = function(X, ths){
  list_id_pep = lapply(unique(X$sequence), function(x){
    id = which(X$sequence == x)
    pep_values = X[id, "PEP"]
    keep = rep(T, length(pep_values))
    
    while (any(keep)) {
      equal = abs(pep_values[keep][1] - pep_values[keep]) < ths
      pep_values[keep][equal] = median(pep_values[keep][equal])
      keep[keep][equal] = F
    }
    list(id, pep_values)
  })
  
  for (i in seq_len(length(list_id_pep))) {
    x = list_id_pep[[i]]
    X[x[[1]], "PEP"] = x[[2]]
  }
  
  COLLAPSED_COUNTS_PEP = aggregate.data.frame(X$Y, by = list(X$PEP, X$sequence), FUN = sum)
  colnames(COLLAPSED_COUNTS_PEP) = c("PEP", "sequence", "Y")
  
  # add EC
  PEPTIDES_EC = X[!duplicated(X$sequence), c("sequence", "EC")]
  COLLAPSED_COUNTS_PEP = merge(COLLAPSED_COUNTS_PEP, PEPTIDES_EC, by = "sequence")
  COLLAPSED_COUNTS_PEP$sequence = NULL
  
  colnames(COLLAPSED_COUNTS_PEP) = c("PEP", "Y", "EC")
  COLLAPSED_COUNTS_PEP
}