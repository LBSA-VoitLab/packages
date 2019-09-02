#' getSignTFs
#'
#' This function gets a list of significant TFs from Trust analysis from specified genes
#' @param gene_set_groups
#' @keywords GSEA gmt
#' @export
#' @examples
#' getSignTFs()

#get Selected TFs from gene sets
getSignTFs <- function(genes){
  set.seed(2235)
  setMode(2)
  trrust_sub <- filterGenes(trrust,genes)
  if(dim(trrust_sub)[1] == 0) return(0)
  nodes_sub <- getNodes(trrust_sub)
  edges_sub <- getEdges (trrust_sub)
  nodes_sub <- getWeights(nodes_sub,edges_sub)
  nodes_sub <- getClusters(nodes_sub,edges_sub)
  #edges_sub <- getPMIDs(edges_sub,trrust_sub)
  #edges_sub <- getAction(edges_sub,trrust_sub)

  pVal <- c()
  for (thisTF in unique(trrust_sub$TF_gene)){
    ids = unique(trrust$nonTF_gene[trrust$TF_gene == thisTF])
    x = sum(ids %in% genes)
    n = length(ids)
    X = sum(unique(trrust$nonTF_gene) %in% genes)
    N = length(unique(trrust$nonTF_gene))

    p <- phyper(x - 1, X, N - X, n, lower.tail = FALSE)
    pVal <- c(pVal,p)
  }
  df <- data.frame(tfs = unique(trrust_sub$TF_gene),
                   pV = pVal,
                   qV = p.adjust(pVal, method = "fdr"))
  return(df[df$pV < 0.05,])
}

#take genes as input and give out a table of genes and their TFs
getGenesSubsetByTf <- function(genes){
#  print(genes)
  tfs_table = getSignTFs(genes)
#  print(tfs_table)
  tfs =tfs_table$tfs
#  print(tfs)
  setMode(2)
  trrust_sub <- filterGenes(trrust,genes)
  setMode(1)
  trrust_sub.2 = list()
  for (tf in tfs){
    this_trrust = filterGenes(trrust,tf)
    trrust_sub.2=rbind(trrust_sub.2,this_trrust[,c(1,2)])
  }
  return(trrust_sub.2)
}
