load_tpm = function(protein_name, mRNAfile){
	tpm = fread(file = mRNAfile, sep = '\t', header = TRUE)
	matches = match(protein_name, tpm$isoname)
	tpm = tpm$tpm[matches]
	
	if(sum(is.na(tpm)) > 0){
	  not_matched = which(is.na(tpm))
	  print("!!! WARNING !!!")
	  print(glue("Found {length(matches) - length(not_matched)} correspondances between transcripts and protein isoforms.
	              The total number of protein isoforms is {length(protein_name)}. For {length(not_matched)} trascript(s) we set the tpm equal to 0.")
	        )
	  print("If the number of correspondances is relatively high,
	        check if the format name of the protein isoforms is consistent with those of the transcripts.")
	  tpm[is.na(tpm)] = 0
	}
	
	# normalize TPMs:
	tpm = tpm/sum(tpm) * 10^6
	
	tpm
}
