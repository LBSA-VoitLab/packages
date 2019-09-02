#' diffExp
#'
#' This function allows to perform deseq2 diff exp analysis 
#' @param gene_set_groups 
#' @keywords GSEA gmt
#' @export
#' @examples
#' diffExp()


library("DESeq2")
library("tximport")
library("readr")
library("tximportData")

diffExp <- function(df,samples,des,cont){
  a <- df[!rowSums(is.na(df)) > 0,]
  dds_data <- DESeqDataSetFromMatrix(countData = a,
                                     colData = samples,
                                     design = ~Exp)
  design(dds_data) = des
  mm=model.matrix(des,colData(dds_data))
  mm=mm[,colSums(mm)>0]
  dds=DESeq(dds_data,full = mm,betaPrior = F)
  
  
  return(list("dds"=dds,"res"=res))
}