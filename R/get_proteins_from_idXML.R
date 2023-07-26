get_proteins_from_idXML = function(file){
  #PG = read_xml(file)
  idXML = data.table::fread(file, sep = NULL, header = FALSE)
  keep = rep(FALSE, nrow(idXML))
  keep = keep + (substr(idXML$V1, 1, 14) == "\t\t\t<ProteinHit")
  idXML = idXML[keep > 0, ]
  
  #PROTEINS = get_records(PG, "//ProteinHit")
  
  #data.frame(id = get_attributes(PROTEINS, "id=\"", "\" accession"),
   #          isoform = get_attributes(PROTEINS, "accession=\"", "\" score="))
  data.frame(id = get_attributes(idXML$V1, "id=\"", "\" accession"),
             isoform = get_attributes(idXML$V1, "accession=\"", "\" score="))
}