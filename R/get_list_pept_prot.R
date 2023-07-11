get_list_pept_prot = function(pept_unique_df, groups, N_peptides_per_protein, pp, params){
  pept_unique_not_in_prot = !(pept_unique_df$EC_numeric %in% unique(unlist(lapply(groups, function(x){x$proteins}))))
  pept_prot = unique(pept_unique_df[pept_unique_not_in_prot,]$EC_numeric)
  pept_prot_df = pept_unique_df[pept_unique_not_in_prot,]
  
  list_pept_prot = lapply(pept_prot, function(x){
    sub_df = pept_prot_df[pept_prot_df$EC_numeric == x, ]
    list(EC_numeric = as.list(1), Y = 0, PEP = 1, M = 1,
         EC_numeric_unique = rep(1, nrow(sub_df)),
         Y_unique = sub_df$Y,
         PEP_unique = sub_df$PEP,
         M_unique = nrow(sub_df),
         N_peptides_per_protein = N_peptides_per_protein[x],
         pp = pp[x], N = 1, K = params$K,
         burn_in = params$burn_in,
         thin = params$thin)
  })
  list(list_pept_prot, pept_prot)
}