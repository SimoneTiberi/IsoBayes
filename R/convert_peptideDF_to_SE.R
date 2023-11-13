convert_peptideDF_to_SE = function(peptideDF, input_type, PEP, FDR_thd,
                                   protein_name, id_openMS = NULL){
  
  Y = matrix(peptideDF$Y)
  sel_col = colnames(peptideDF)[-which(colnames(peptideDF) == "Y")]
  feature_data = DataFrame(peptideDF[, sel_col])
  colnames(feature_data) = sel_col
  
  SE = SummarizedExperiment(assays=list(Y=t(Y)), 
                            colData=feature_data,
                            metadata = list(PEP = PEP,
                                            FDR_thd = FDR_thd,
                                            input_type = input_type,
                                            protein_name = protein_name,
                                            id_openMS = id_openMS))
  SE
}