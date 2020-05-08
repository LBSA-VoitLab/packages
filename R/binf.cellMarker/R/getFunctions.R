#' getCellNames
#'
#' This function gets values from marker genes
#' @param gene_set_groups
#' @keywords GSEA gmt
#' @export
#' @examples
#' getCellNames()


#data(cell_markers_db)
#get cell name from the db
getCellNames <- function(gene){
  marker_filter = sapply(all_marker_genes,function(x) ifelse(sum(unlist(x) %in% gene) > 0, TRUE, FALSE))

  return(levels(factor(cell_markers_db[marker_filter,"cellName"])))
}

#get all gene markers
getAllGeneMarkers <- function(){
  return(unique(sapply( unlist(all_marker_genes), paste0, collapse="")))
}

#get all gene markers array with specified db
getAllMarkersArray <- function(db){
  markers_db = getAllMarkersList(db)
  return(unique(sapply( unlist(markers_db), paste0, collapse="")))
}

#get all gene markers with specified db
getAllMarkersList <- function(db){
  all_marker_genes = db$geneSymbol
  markers_db = lapply(all_marker_genes,function(x) unlist(lapply(strsplit(levels(factor(x)),split = ","),trimws)))
  return(markers_db)
}

#get sub array
getSubDb <- function(genes){
  marker_filter = sapply(all_marker_genes,function(x) ifelse(sum(unlist(x) %in% genes) > 0, TRUE, FALSE))
  return(cell_markers_db[marker_filter,])
}

#get filtered array
getFiltDb <- function(genes,category,value){
  sub_db = getSubDb(genes)
  return(sub_db[sub_db[,category] == value,])
}

#get sub db by category
getSubDbCategory <- function(category,value){
  return(cell_markers_db[cell_markers_db[,category] %in% value,])
}

#get sub db by category 2
getSubDbCategory.2 <- function(db,category,value){
  return(db[db[,category] %in% value,])
}

#getCell type from db
getCellType <- function(db){
  return(as.character(db$cellName))
}

#get enrichment of of cellnames by tissue type
marker_enrich = function (gene_list,category,value,category2,value2)
{
  db = getSubDbCategory(category,value)
  db = getSubDbCategory.2(db,category2,value2)
  all_markers_list = getAllMarkersList(db)
  all_markers_array= getAllMarkersArray(db)
  cell_type=getCellType(db)
  #print(gene_list)

  gene_list=intersect(gene_list,all_markers_array)

  #print(gene_list)

  pop.size=length(all_markers_array)
  samp.size=length(gene_list)

  p.val=rep(0,length(all_markers_list))
  samp.hits=rep(0,length(all_markers_list))
  pop.hits=rep(0,length(all_markers_list))
  genes=list()

  for (k in 1:length(all_markers_list))
  {
    samp.hits[k]=length(intersect(gene_list,all_markers_list[[k]]))
    pop.hits[k]=length(all_markers_list[[k]])
    p.val[k]=phyper(samp.hits[k], pop.hits[k], pop.size-pop.hits[k], samp.size,lower.tail = FALSE)+dhyper(samp.hits[k], pop.hits[k], pop.size-pop.hits[k], samp.size)
    genes[[k]]=intersect(gene_list,all_markers_list[[k]])
  }

  ix=order(p.val)

  enrich=data.frame(Cell.type=cell_type, p.value=p.val,overlap=samp.hits,signature=pop.hits)
  enrich=enrich[ix,]

  result=list(enrichments=enrich,genes=genes[ix])

}
