get_records = function(pg, record){
  as.character(xml_find_all(pg, record))
}

get_attributes = function(record, begins, ends, isNumeric = FALSE, multiple = FALSE){
  out = sapply(record, function(x){gsub(paste0(".*", begins), "", x)})
  if(isNumeric){
    out = sapply(out, function(x){as.numeric(gsub(paste0(ends, ".*"), "", x))})
  }else{
    out = sapply(out, function(x){gsub(paste0(ends, ".*"), "", x)})
  }
  names(out) = NULL
  
  out
}