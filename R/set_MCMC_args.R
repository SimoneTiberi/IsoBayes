set_MCMC_args = function(pept_df, pept_unique_df, prot_df, PEP, prior) {
  arguments = list()
  arguments$pept_df = pept_df
  if (PEP) {
    arguments$pept_unique_df = pept_unique_df
    arguments$prot_df = prot_df
    arguments$lib_size = sum(arguments$pept_df$Y * (1 - arguments$pept_df$PEP)) + sum(arguments$pept_unique_df$Y * (1 - arguments$pept_unique_df$PEP))
    arguments$M_unique = nrow(arguments$pept_unique_df) # number of filtered detected unique peptides
    arguments$prior = prior
    arguments$pept_unique_df$EC_numeric = unlist(arguments$pept_unique_df$EC_numeric)
    arguments$pept_unique_df$Y = as.integer(arguments$pept_unique_df$Y)

    if (any(class(arguments$pept_unique_df$EC_numeric) != "integer", class(arguments$pept_unique_df$Y) != "integer")) {
      stop("Error: classes types don't respect C++ classes")
    }
  } else {
    arguments$prot_df = prot_df
    arguments$lib_size = sum(arguments$prot_df$Y_unique) + sum(arguments$pept_df$Y)
  }
  arguments$pept_df$Y = as.numeric(arguments$pept_df$Y)

  if (any(class(arguments$pept_df$EC_numeric) != "list", class(arguments$pept_df$Y) != "numeric", class(arguments$prot_df$Y_unique) != "numeric")) {
    stop("Error: classes types don't respect C++ classes")
  }

  arguments$N = nrow(arguments$prot_df) # number of proteins
  arguments$M = nrow(arguments$pept_df) # number of filtered detected peptides
  arguments$prior = prior
  arguments$protein_length = arguments$prot_df$protein_length

  arguments
}
