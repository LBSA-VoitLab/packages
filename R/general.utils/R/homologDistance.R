#' vlookup
#'
#' some function for homolog distance
#' @param this df key value
#' @keywords 
#' @export
#' @examples
#' vlookup(this, df, key, value)


getAverageHomologDistanceForGeneset <- function(geneSetNames,geneDistanceTable,gsea_res_folder){
  gene_set_groups = paste0(gsea_res_folder,"/edb/gene_sets.gmt")
  gsg <- get_geneset_table(gene_set_groups)
  thisSub <-gsg[geneSetNames]
  geneSetDistList <- NULL
  for (i in 1:length(thisSub)){
    #    print(thisSub[[i]])
    thisSetDist <- getAverageHomologDistance(thisSub[[i]],geneDistanceTable)
    #    geneSetDistList <- c(geneSetDistList,thisSetDist)
    geneSetDistList = rbind(geneSetDistList, data.frame(geneSetNames[i],thisSetDist))
    #    print(thisSetDist)
    
  }
  return(geneSetDistList)
  #  return(list("gene" <- geneSetNames,"dist" <- geneSetDistList))
}

getAverageHomologDistance <- function(geneSet,geneDistanceTable){
  ind <- geneDistanceTable$gene %in% geneSet
  distList <- geneDistanceTable$dist[ind]
  #  print(geneDistanceTable$gene[ind])
  #  print(geneDistanceTable$dist[ind])
  #  print(distList)
  return(mean(distList))
}

convertHomologTable2DistTable <- function(all_homologs.filter){
  geneDistanceTable <- c()
  geneDistanceTable[["gene"]] <- as.character(all_homologs.filter$gene.name.Mm)
  geneDistanceTable[["dist"]] <- all_homologs.filter$dist
  return(geneDistanceTable)
}