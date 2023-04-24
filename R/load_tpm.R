load_tpm = function(protein_name, mRNAfile){
	tpm = fread(file = mRNAfile, sep = '\t', header = TRUE)
	matches = match(protein_name, tpm$isoname)
	tpm = tpm$tpm[matches]
	tpm[is.na(tpm)] = 0

	# normalize TPMs:
	tpm = tpm/sum(tpm) * 10^6
	
	tpm
}
