aggregate_components_pep = function(components, ncores){
  nc = length(components)
  for(i in 1:ncores){
    for (j in 1:floor(nc/ncores)) {
      if(i < length(components)){
        a = lapply(components[[i+1]]$EC_numeric, function(x){x+components[[i]]$N})
        components[[i]]$EC_numeric = append(components[[i]]$EC_numeric, a)
        components[[i]]$Y = c(components[[i]]$Y, components[[i+1]]$Y)
        components[[i]]$PEP = c(components[[i]]$PEP, components[[i+1]]$PEP)
        components[[i]]$M = components[[i]]$M + components[[i+1]]$M
        components[[i]]$EC_numeric_unique = c(components[[i]]$EC_numeric_unique,
                                              components[[i+1]]$EC_numeric_unique + components[[i]]$N)
        components[[i]]$Y_unique = c(components[[i]]$Y_unique, components[[i+1]]$Y_unique)
        components[[i]]$PEP_unique = c(components[[i]]$PEP_unique, components[[i+1]]$PEP_unique)
        components[[i]]$M_unique = components[[i]]$M_unique + components[[i+1]]$M_unique
        components[[i]]$protein_length = c(components[[i]]$protein_length, components[[i+1]]$protein_length)
        components[[i]]$pp = c(components[[i]]$pp, components[[i+1]]$pp)
        components[[i]]$N = components[[i]]$N + components[[i+1]]$N
        components[[i+1]] = NULL
      }
    }
  }
  components
}