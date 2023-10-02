unique_peptides = function(EC_numeric) {
    EC_length = vapply(EC_numeric, length, FUN.VALUE = numeric(1))

    sel_unique = EC_length == 1 # peptides with 1 protein only
    message(glue("Percentage of unique peptides:
                 {round(mean(sel_unique), 4)*100}%"))

    sel_unique
}