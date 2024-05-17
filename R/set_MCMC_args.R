set_MCMC_args = function(pept_df, pept_unique_df, prot_df, PEP, prior) {
    arguments = list()
    arguments$pept_df = pept_df
    if (PEP) {
        arguments$pept_unique_df = pept_unique_df
        arguments$prot_df = prot_df
        sum_pept = sum(arguments$pept_df$Y * (1 - arguments$pept_df$PEP))
        sum_pept_unique = sum(arguments$pept_unique_df$Y * (1 - arguments$pept_unique_df$PEP))
        arguments$lib_size = sum_pept + sum_pept_unique
        arguments$M_unique = nrow(arguments$pept_unique_df)
        arguments$pept_unique_df$EC_numeric = unlist(arguments$pept_unique_df$EC_numeric)
        arguments$pept_unique_df$Y = as.integer(arguments$pept_unique_df$Y)

        if (any(!is.integer(arguments$pept_unique_df$EC_numeric),
                !is.integer(arguments$pept_unique_df$Y))) {
            stop("classes types don't respect C++ classes")
        }
    } else {
        arguments$prot_df = prot_df
        arguments$lib_size = sum(arguments$prot_df$Y_unique) + sum(arguments$pept_df$Y)
    }
    arguments$pept_df$Y = as.numeric(arguments$pept_df$Y)

    if (any(!is.list(arguments$pept_df$EC_numeric),
            !is.numeric(arguments$pept_df$Y),
            !is.numeric(arguments$prot_df$Y_unique))) {
        stop("classes types don't respect C++ classes")
    }

    arguments$N = nrow(arguments$prot_df) # number of proteins
    arguments$M = nrow(arguments$pept_df) # number of filtered detected peptides
    arguments$prior = prior
    arguments$protein_length = arguments$prot_df$protein_length

    arguments
}