get_overall_abundance = function(Y_shared, Y_unique, prot_unique) {
    if (is.null(Y_unique)) {
        sum(Y_shared$Y) + sum(prot_unique)
    } else {
        sum(Y_shared$Y * (1 - Y_shared$PEP)) +
            sum(Y_unique$Y * (1 - Y_unique$PEP))
    }
}