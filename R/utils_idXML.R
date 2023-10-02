get_attributes = function(record, begins, ends, isNumeric = FALSE) {
    out = vapply(record, function(x) {
        gsub(paste0(".*", begins), "", x)
    }, FUN.VALUE = character(1))
    if (isNumeric) {
        out = vapply(out, function(x) {
            as.numeric(gsub(paste0(ends, ".*"), "", x))
        }, FUN.VALUE = numeric(1))
    } else {
        out = vapply(out, function(x) {
            gsub(paste0(ends, ".*"), "", x)
        }, FUN.VALUE = character(1))
    }
    names(out) = NULL

    out
}