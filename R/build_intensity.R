build_intensity = function(z, path_to_peptides_psm, variables) {
    sel_intensities = grep("Intensity_", colnames(z))
    # ADD the signal across the k intensities:
    intensity = rowSums(z[, sel_intensities])
    message(glue("Peptides with intensity equal to 0:
                 {round(mean(intensity == 0), 5)*100}%"))

    # AllQuantifiedPeptides.tsv does not have information about ECs
    # We load the "PSM" file AllPeptides.psmtsv below and match intensities according to "Base Sequence"
    # "Base Sequence" is the sequence of the peptide
    x = as.data.frame(fread(path_to_peptides_psm, select = variables))
    matches = match(x$`Base Sequence`, z$`Base Sequence`)
    message(glue("NA's intensity: {round(mean(is.na(matches)), 5)*100}%"))

    INTENSITY = intensity[matches]
    INTENSITY[is.na(INTENSITY)] = 0

    list(y = INTENSITY, x = x)
}