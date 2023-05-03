#' @useDynLib SIMBA, .registration = TRUE
#' @importFrom methods formalArgs
#' @importFrom stats aggregate.data.frame
#' @importFrom Biostrings fasta.seqlengths
#' @importFrom data.table fread
#' @importFrom data.table setnames
#' @importFrom glue glue
#' @importFrom Rcpp sourceCpp
#' @importFrom xml2 read_xml
#' @importFrom xml2 xml_find_all
#' @importFrom DescTools Mode
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG '%dorng%'
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom foreach foreach
#' @importFrom iterators iter
#' @importFrom stats median
#' @importFrom stats na.omit
#' @importFrom stats quantile
NULL
