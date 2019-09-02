
#' createTableAndSelect
#'
#' This function allows to perform deseq2 diff exp analysis 
#' @param gene_set_groups 
#' @keywords GSEA gmt
#' @export
#' @examples
#' createTableAndSelect()



#from deseq results to table and selected table
createTableAndSelect <- function(res,p = NULL,fc = NULL){
  if(is.null(p)) p <- 0.01
  if(is.null(fc)) fc <- 1
  
  tab.1 <- data.frame(fc = res$log2FoldChange,pVal = res$pvalue,row.names = row.names(res))
  tab <- tab.1[!is.na(tab.1$pVal),]
  sel <- tab[tab$pVal < p & abs(tab$fc) > fc,]
  return(list(tab = tab,sel = sel))
}
