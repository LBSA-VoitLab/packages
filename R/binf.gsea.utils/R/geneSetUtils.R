#' get_geneset_table
#'
#' This function allows to create a list of list for GSEA geneset (.gmt) file
#' @param gene_set_groups 
#' @keywords GSEA gmt
#' @export
#' @examples
#' get_geneset_table()
 
 
#gene_set_groups = "Z:/projects-bmed/projects/AHR/E6v7_hallmark/E6v7_hallmark.GseaPreranked.1553106928010/edb/gene_sets.gmt" 

get_geneset_table <- function(gene_set_groups){
  gsg <- readLines(gene_set_groups)
  df <- list()
  for(i in c(1:length(gsg))){
    line <- unlist(strsplit(gsg[i],split = "\t"))
    line = line[line != ""]
    df[[toupper(as.character(line[1]))]] <-line[3:length(line)]
  }
  return(df)
}

#from gsea result folder
get_geneset_table.2 <- function(gsea_folder){
  gene_set_groups = paste0(gsea_folder,"/edb/gene_sets.gmt")
  return(get_geneset_table(gene_set_groups))
}

#subset of geneset_table by geneset names
subset_geneset_table <- function(geneset_table, geneset_names){
  sub_geneset_table = geneset_table[toupper(names(gsg)) %in% geneset_names]
  return(sub_geneset_table)
}

#from gsea results folder get top n pos geneset names and genesets
get_pos_genesets <- function(gsea_res_folder,n = NULL){
  if(is.null(n)) n = 10
  posData <- get_all_pos_list(gsea_res_folder)
  sortedPosData = posData[order(-posData$NES),]
  top_n_sortedPosData = sortedPosData[1:n,]
  return(top_n_sortedPosData)
}

#from gsea results folder get all pos geneset names
get_all_pos_list <- function(gsea_res_folder){
  posList = list.files(path=gsea_res_folder,pattern="gsea_report_for_.*_pos_.*.xls")
  posData <- read.csv2(paste0(gsea_res_folder,"/",posList), header=TRUE,sep="\t")
  return(posData)
}

#from gsea results folder get top n neg geneset names and genesets
get_neg_genesets <- function(gsea_res_folder,n = NULL){
  if(is.null(n)) n = 10
  negData <- get_all_neg_list(gsea_res_folder)
  sortedNegData = negData[order(negData$NES),]
  top_n_sortedNegData = sortedNegData[1:n,]
  return(top_n_sortedNegData)
}

#from gsea results folder get all neg geneset names
get_all_neg_list <- function(gsea_res_folder){
  negList <- list.files(path=gsea_res_folder,pattern="gsea_report_for_.*_neg_.*.xls")
  negData <- read.csv2(paste0(gsea_res_folder,"/",negList), header=TRUE,sep="\t")
  return(negData)
}
