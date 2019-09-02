library(	data.table)


pre <- function(){
  trrust_old <- fread("Z:/projects-bmed/projects/packages/R/binf.trrust/data/trrust_rawdata.txt", header=FALSE, stringsAsFactors = FALSE, showProgress = FALSE, data.table = FALSE)
  trrust <- fread("Z:/projects-bmed/projects/packages/R/binf.trrust/data/trrust_rawdata.human.tsv", header=FALSE, stringsAsFactors = FALSE, showProgress = FALSE, data.table = FALSE)

  colnames(trrust)[1] <- "TF_gene"
  colnames(trrust)[2] <- "nonTF_gene"
  colnames(trrust)[3] <- "interaction"
  colnames(trrust)[4] <- "PMID"

  save(trrust, file = "Z:/projects-bmed/projects/packages/R/binf.trrust/data/trrust.RData")
}
#Search TFs by known genes
get_tf_by_genes <- function(){
  edges <- getEdges (trrust)
  nodes <- getNodes(trrust)
}

#Search genes by known Tfs
get_genes_by_tf <- function(){

}

###
#Get Edges
###

#x is a dataframe where you imported your TRRUST network into

getEdges <-function(x){
  igraph_object <- graph.data.frame(x)
  edges <- data.frame(get.edgelist(igraph_object))

  #Fix the Column Names
  colnames(edges)[1] <- "from"
  colnames(edges)[2] <- "to"


  return(edges)

}
