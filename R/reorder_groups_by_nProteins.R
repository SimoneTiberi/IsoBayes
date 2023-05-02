reorder_groups_by_nProteins = function(groups, params){
  groups_weight = unlist(lapply(groups, function(x){length(x$proteins)}))
  new_order = sort(groups_weight, index.return = TRUE, decreasing = TRUE)$ix
  distributed_index = c()
  for (inc in seq_len(params$n_cores)) {
    distributed_index = c(distributed_index, seq(0, length(new_order), params$n_cores) + inc)
  }
  new_order = na.omit(new_order[distributed_index])
  
  groups[new_order]
}
