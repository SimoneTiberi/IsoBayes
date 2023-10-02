unique_protein_abundance = function(y, EC, protein_name) {
    EC_numeric = convert_EC_to_num(strsplit(EC, split = "\\|"), protein_name)
    sel_unique = unique_peptides(EC_numeric)

    Y_unique = rep(0, length(protein_name))
    for (i in which(sel_unique)) {
        Y_unique[unlist(EC_numeric[i])] = Y_unique[unlist(EC_numeric[i])] + y[i]
    }

    list(Y_unique = Y_unique, sel_unique = sel_unique, EC_numeric = EC_numeric)
}