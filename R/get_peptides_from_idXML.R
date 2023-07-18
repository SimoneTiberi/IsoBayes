get_peptides_from_idXML = function(file, pep, FDR_thd){
  PG = read_xml(file)
  PEPTIDES = get_records(PG, "//PeptideHit")

  if(FDR_thd < 1){
    FDR = get_attributes(PEPTIDES, "q-value\" value=\"", "\"/>\n</PeptideHit>", 
                         isNumeric = TRUE)
    if(mean(is.na(FDR)) > 0){
      stop(glue("{mean(is.na(FDR)*100}% of FDRs not found. Check the .idXML file: for each 'PeptideHit', 'UserParam' q-value should be present.
                For more details on how to generate an .idXML, see the vignettes."))
    }
  }else{
    FDR = NULL
  }
  if(pep){
    PEP = get_attributes(PEPTIDES, "score=\"", "\" sequence", isNumeric = TRUE)
  }else{
    PEP = NULL
  }
  
  PEPTIDES = data.frame(Y = 1,
                        EC = get_attributes(PEPTIDES, "protein_refs=\"", "\">\n  <UserParam"),
                        target_decoy = get_attributes(PEPTIDES, "name=\"target_decoy\" value=\"",
                                                      "\"/>\n  <UserParam type"), 
                        sequence = get_attributes(PEPTIDES, "sequence=\"", "\" charge=")
                        )
  PEPTIDES$PEP = PEP
  PEPTIDES$FDR = FDR  
  rm(PG)
  
  PEPTIDES = PEPTIDES[PEPTIDES$target_decoy == "target", ]
  if(FDR_thd < 1){
    PEPTIDES = PEPTIDES[PEPTIDES$FDR < FDR_thd, ]
  }
  PEPTIDES$target_decoy = NULL
  PEPTIDES$FDR = NULL
  
  if(pep){
    # set equal pep if difference < 1% and sequence is equal
    equal_pep = function(x){
      PEPTIDES[x[[1]], "PEP"] = x[[2]]
    }
    
    temp_PEPTIDES = PEPTIDES
    list_id_pep = lapply(unique(temp_PEPTIDES$sequence), function(x){
      id = which(temp_PEPTIDES$sequence == x)
      pep_values = temp_PEPTIDES[id, "PEP"]
      keep = rep(TRUE, length(pep_values))
      
      while (any(keep)) {
        equal = abs(pep_values[keep][1] - pep_values[keep]) < 0.01
        pep_values[keep][equal] = median(pep_values[keep][equal])
        keep[keep][equal] = FALSE
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
