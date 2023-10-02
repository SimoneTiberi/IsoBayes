collapse_pept_w_equal_EC = function(pept_df, pep) {
    if (pep) {
        aggregate_sim_pep(pept_df, ths = 0.01)
    } else {
        pept_df$sequence = NULL
        pept_df = aggregate.data.frame(pept_df$Y, by = list(pept_df$EC),
                                       FUN = sum)
        colnames(pept_df) = c("EC", "Y")
        pept_df
    }
}