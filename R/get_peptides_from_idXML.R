get_peptides_from_idXML = function(file, pep){
  PG = read_xml(file)
  PEPTIDES = get_records(PG, "//PeptideHit")
  PEPTIDES = data.frame(Y = 1,
                        EC = get_attributes(PEPTIDES, "protein_refs=\"", "\">\n  <UserParam"),
                        target_decoy = get_attributes(PEPTIDES, "name=\"target_decoy\" value=\"", "\"/>\n  <UserParam type"),
                        PEP = get_attributes(PEPTIDES, "score=\"", "\" sequence", isNumeric = TRUE),
                        sequence = get_attributes(PEPTIDES, "sequence=\"", "\" charge=")
                        )
  rm(PG)
  
  PEPTIDES = PEPTIDES[PEPTIDES$target_decoy == "target", ]
  PEPTIDES$target_decoy = NULL
  
  if(pep){
    # set equal pep if difference < 1% and sequence is equal
    equal_pep = function(x){
      PEPTIDES[x[[1]], "PEP"] = x[[2]]
    }
    
    temp_PEPTIDES = PEPTIDES
    list_id_pep = lapply(unique(temp_PEPTIDES$sequence), function(x){
      id = which(temp_PEPTIDES$sequence == x)
      pep_values = temp_PEPTIDES[id, "PEP"]
      keep = rep(T, length(pep_values))
      
      while (any(keep)) {
        equal = abs(pep_values[keep][1] - pep_values[keep]) < 0.01
        pep_values[keep][equal] = median(pep_values[keep][equal])
        keep[keep][equal] = F
      }
      list(id, pep_values)
    })
    
    for (i in seq_len(length(list_id_pep))) {
      x = list_id_pep[[i]]
      temp_PEPTIDES[x[[1]], "PEP"] = x[[2]]
    }
    
    COLLAPSED_COUNTS = aggregate.data.frame(temp_PEPTIDES$Y, by = list(temp_PEPTIDES$PEP, temp_PEPTIDES$sequence), FUN = sum)
    rm(temp_PEPTIDES, list_id_pep)
    colnames(COLLAPSED_COUNTS) = c("PEP", "sequence", "Y")
    
    PEPTIDES_EC = PEPTIDES[!duplicated(PEPTIDES$sequence), c("sequence", "EC")]
    rm(PEPTIDES)
    COLLAPSED_COUNTS = merge(COLLAPSED_COUNTS, PEPTIDES_EC, by = "sequence")
    rm(PEPTIDES_EC)
    COLLAPSED_COUNTS$EC = gsub(" ", "|", COLLAPSED_COUNTS$EC)
  }else{
    COLLAPSED_COUNTS = aggregate.data.frame(PEPTIDES$Y, by = list(PEPTIDES$EC), FUN = sum)
    colnames(COLLAPSED_COUNTS) = c("EC", "Y")
    COLLAPSED_COUNTS$EC = gsub(" ", "|", COLLAPSED_COUNTS$EC)
    rm(PEPTIDES)
  }
  
  COLLAPSED_COUNTS
}
