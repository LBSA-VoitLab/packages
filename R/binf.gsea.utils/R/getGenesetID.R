#' get_geneset_ID Function
#'
#' This function allows you to search geneset id from name
#' @param geneset_name name of geneset to be searched
#' @keywords barplot gsea
#' @export
#' @examples
#' get_geneset_ID()

get_geneset_ID <- function(geneset_name){
  goterms <- Term(GOTERM)
  c <- ClosestMatch2(geneset_name,goterms)
  return(c)
}


ClosestMatch2 = function(string, stringVector){
  stringVector[amatch(string, stringVector, maxDist=Inf)]
} 