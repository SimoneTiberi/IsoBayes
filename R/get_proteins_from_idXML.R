get_proteins_from_idXML = function(idXML) {
    keep = rep(FALSE, nrow(idXML))
    keep = keep + (substr(idXML$V1, 1, 14) == "\t\t\t<ProteinHit")
    idXML = idXML[keep > 0, ]

    data.frame(
        id = get_attributes(idXML$V1, "id=\"", "\" accession"),
        isoform = get_attributes(idXML$V1, "accession=\"", "\" score=")
    )
}
