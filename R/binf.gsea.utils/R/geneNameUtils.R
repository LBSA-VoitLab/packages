#' getMapForList
#'
#' This function allows to create a mapping for 2 sets of gene names based on NCBI gene_info file
#' @param gene_set_groups 
#' @keywords GSEA gmt
#' @export
#' @examples
#' getMapForList()


#gene_set_groups = "Z:/projects-bmed/projects/AHR/E6v7_hallmark/E6v7_hallmark.GseaPreranked.1553106928010/edb/gene_sets.gmt" 
library(foreach)
library(doParallel)
library(rlist)

getMapForList <- function(lst1,lst2,gene_info){
  no_cores <- detectCores() - 1
  cl<-makeCluster(no_cores)
  registerDoParallel(cl)
  
  sym <- gene_info$Symbol
  syn <- strsplit(as.character(gene_info$Synonyms),"\\|")
  id <- paste0("LOC",as.character(gene_info$GeneID))
  tab <- data.frame()
  #  tab <- foreach(i = 1:length(lst1), .combine = rbind, .export = "getMapForGene") %dopar% getMapForGene(lst1[i],lst2,sym,syn,id)
  tab <- foreach(i = 1:length(lst1), .combine = rbind, .export = c("getMapForGene","list.is")) %dopar% getMapForGene(lst1[i],lst2,sym,syn,id)
  #  names(tab) <- c("E06","E07")
  return(tab)
}
getMapForGene <- function(gn,lst2,sym,syn,id){
  #  candidate_subset <- (sym %in% gn) | grepl(paste0("^",gn,"$"),syn) | (id %in% gn)
  candidate_subset <- (sym %in% gn) | list.is(syn,any(. == gn)) | (id %in% gn)
  candidates <- unique(c(as.character(sym[candidate_subset]), unlist(syn[candidate_subset]),as.character(id[candidate_subset])))
  selected <- candidates[candidates %in% lst2]
  if (length(selected) > 0){
    print(cbind(gn,selected))
    return(cbind(gn,selected))
  }
  else{
    #      print(c("N",gn))
    return(c(gn,"0"))
  }
}