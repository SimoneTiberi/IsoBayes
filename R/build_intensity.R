build_intensity = function(z, path_to_peptides_psm){
  sel_intensities = grep("Intensity_", colnames(z))
  # ADD the signal across the 6 intensities:
  intensity = rowSums(z[, sel_intensities])
  print("Peptides with intensity == 0:") ; print(mean(intensity == 0))
  
  # AllQuantifiedPeptides.tsv does not have information about ECs
  # We load the "PSM" file AllPeptides.psmtsv below and match intensities according to "Base Sequence"
  # "Base Sequence" is the sequence of the peptide
  x = as.data.frame(fread(path_to_peptides_psm))
  matches = match(x$`Base Sequence`, z$`Base Sequence`)
  print("Na's intensity:") ; print(mean(is.na(matches)))
  
  INTENSITY = intensity[matches]
  INTENSITY[is.na(INTENSITY)] = 0
  
  # numerical ISSUE with intensities: intensities have VERY large values, so sometimes multi-mapping peptides
  # will distribute at least 1 "intensity" to all proteins -> all with have a high probability of being present!
  # 2 solutions:
  # 1) Change threshold to Pr( X > min_count_like_100 )
  # 2) decrease intensities to a "reasonable" scale.
  # below, we follow 2) -> we compute TPM-like intensities (sum of all intensities = 10^5).
  # I choose 10^5 because it's similar to the sum of PSMs in some example datasets we looked at.
  INTENSITY = INTENSITY/sum(INTENSITY) * 10^5
  print("Sum of the intensities") ; print(sum(INTENSITY))
  
  # round intensities to closest integer, BUT we add 0.5 so that very small intensities (between 0 and 0.5) are rounded to 1.
  INTENSITY = round(INTENSITY + 0.5)
  
  list(y = INTENSITY, x = x)
}
