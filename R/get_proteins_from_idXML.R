get_proteins_from_idXML = function(file){
  PG = read_xml(file)
  PROTEINS = get_records(PG, "//ProteinHit")
  
  data.frame(id = get_attributes(PROTEINS, "id=\"", "\" accession"),
             isoform = get_attributes(PROTEINS, "accession=\"", "\" score="))
}
