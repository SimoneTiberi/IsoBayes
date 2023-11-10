load_tpm = function(protein_name, tpm) {
    
    tpm = check_and_load_tpm(tpm)
    matches = match(protein_name, tpm$isoname)
    tpm = tpm$tpm[matches]
    
    if (sum(is.na(tpm)) > 0) {
        not_matched = which(is.na(tpm))
        
        message(glue("Found {length(matches) - length(not_matched)}
        correspondances between transcripts and protein isoforms.
        The total number of protein isoforms is {length(protein_name)}.
        For {length(not_matched)} trascript(s) we set the tpm equal to 0."))
        
        message("If the number of correspondances is relatively high,
                check if the format name of the protein isoforms
                is consistent with those of the transcripts.")
        tpm[is.na(tpm)] = 0
    }
    tpm
}